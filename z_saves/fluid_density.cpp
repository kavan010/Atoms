#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <vector>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
using namespace glm;
using namespace std;

const float c = 299792458.0f / 100000000.0f;    // speed of light in m/s
const float k = 8.9875517923e9f; // Coulomb's constant

// ================= Engine ================= //
struct Engine {
    GLFWwindow* window;
    int WIDTH = 800;
    int HEIGHT = 600;

    Engine () {
        if(!glfwInit()) {
            cerr << "failed to init glfw, PANIC!";
            exit(EXIT_FAILURE);
        };
        window = glfwCreateWindow(WIDTH, HEIGHT, "Atom Sim", nullptr, nullptr);
        if(!window) {
            cerr << "failed to create window, PANIC!";
            glfwTerminate();
            exit(EXIT_FAILURE);
        }
        glfwMakeContextCurrent(window);
        glViewport(0, 0, WIDTH, HEIGHT);
    };
    void run() {
        glClear(GL_COLOR_BUFFER_BIT);
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        double left   = -WIDTH;
        double right  = WIDTH;
        double bottom = -HEIGHT;
        double top    = HEIGHT;
        glOrtho(left, right, bottom, top, -1.0, 1.0);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
    }
    
    void drawCircle(float x, float y, float r, int segments = 100) {
        glBegin(GL_LINE_LOOP);
        for (int i = 0; i < segments; i++) {
            float angle = 2.0f * M_PI * i / segments;
            float dx = r * cosf(angle);
            float dy = r * sinf(angle);
            glVertex2f(x + dx, y + dy);
        }
        glEnd();
    }
    void drawFilledCircle(float x, float y, float r, int segments = 50) {
        glBegin(GL_TRIANGLE_FAN);
        glVertex2f(x, y); // center
        for (int i = 0; i <= segments; i++) {
            float angle = 2.0f * M_PI * i / segments;
            float dx = r * cosf(angle);
            float dy = r * sinf(angle);
            glVertex2f(x + dx, y + dy);
        }
        glEnd();
    }
};
Engine engine;

const int   GRID_W = engine.WIDTH;
const int   GRID_H = engine.HEIGHT;
const float CELL_W = 1.0f;
const float CELL_H = 1.0f;

vector<float> density(GRID_W * GRID_H, 0.0f);
float maxDensity = 1.0f;

// After the Engine struct and global constants...
vec3 densityToColour(float d, float maxD) {
    if (maxD <= 0.0f) maxD = 1.0f;
    float t = glm::clamp(d / maxD, 0.0f, 1.0f);
    float r, g, b;

    if (t < 0.33f) {
        // black → red
        float k = t / 0.33f;
        r = k;
        g = 0.0f;
        b = 0.0f;
    } else if (t < 0.66f) {
        // red → yellow
        float k = (t - 0.33f) / 0.33f;
        r = 1.0f;
        g = k;
        b = 0.0f;
    } else {
        // yellow → white
        float k = (t - 0.66f) / 0.34f;
        r = 1.0f;
        g = 1.0f;
        b = k;
    }

    return vec3(r, g, b);
}
void calculateDensityMap(const vector<vec2>& particles) {
    // 1. Reset the density map and max density
    fill(density.begin(), density.end(), 0.0f);
    maxDensity = 0.0f;

    // Define the influence radius (e.g., 25 units in your coordinate system)
    const float INFLUENCE_RADIUS = 200.0f;
    const float INFLUENCE_RADIUS_SQ = INFLUENCE_RADIUS * INFLUENCE_RADIUS;
    
    // Scale factor to make the density values manageable/visible
    const float DENSITY_SCALE = 1.0f; 

    // 2. Iterate through all particles and contribute to the grid
    for (const vec2& p : particles) {
        // Map particle position to grid coordinates (assuming center is 0,0)
        // Convert from OpenGL coordinates (-WIDTH to WIDTH) to Grid indices (0 to GRID_W-1)
        int center_i = (int)round(p.x + GRID_W) / 2; // (p.x + 800) / 2 -> 0 to 800
        int center_j = (int)round(p.y + GRID_H) / 2; // (p.y + 600) / 2 -> 0 to 600

        // Determine the bounds of the grid cells to update (a box around the particle)
        int min_i = std::max(0, center_i - (int)ceil(INFLUENCE_RADIUS));
        int max_i = std::min(GRID_W - 1, center_i + (int)ceil(INFLUENCE_RADIUS));
        int min_j = std::max(0, center_j - (int)ceil(INFLUENCE_RADIUS));
        int max_j = std::min(GRID_H - 1, center_j + (int)ceil(INFLUENCE_RADIUS));

        // Iterate over the nearby grid cells
        for (int j = min_j; j <= max_j; ++j) {
            for (int i = min_i; i <= max_i; ++i) {
                // Get the center position of the grid cell
                // Convert grid index back to OpenGL coordinates
                float cell_x = i * 2.0f - GRID_W;
                float cell_y = j * 2.0f - GRID_H;

                // Calculate the squared distance from particle center (p) to cell center (cell_x, cell_y)
                float dx = p.x - cell_x;
                float dy = p.y - cell_y;
                float distSq = dx * dx + dy * dy;

                if (distSq < INFLUENCE_RADIUS_SQ) {
                    // Use a simple quadratic kernel (smooth falloff)
                    // You could use a Gaussian kernel for smoother results, but this is simpler
                    float dist = sqrt(distSq);
                    float influence = (INFLUENCE_RADIUS - dist) / INFLUENCE_RADIUS;
                    // Contribution = (1 - (dist / radius))^2 
                    float kernel_weight = influence * influence * DENSITY_SCALE; 
                    
                    int index = j * GRID_W + i;
                    density[index] += kernel_weight;
                    maxDensity = std::max(maxDensity, density[index]);
                }
            }
        }
    }
}
// Function to draw the density map as colored squares
void drawDensityMap() {
    // 1. Set up for drawing quads
    glBegin(GL_QUADS);

    // 2. Iterate over the entire grid
    for (int j = 0; j < GRID_H; ++j) {
        for (int i = 0; i < GRID_W; ++i) {
            int index = j * GRID_W + i;
            float d = density[index];
            vec3 color = densityToColour(d, maxDensity);
            
            // Set the color for the current cell (from blue for low to red for high)
            glColor3f(color.r, color.g, color.b);

            // Calculate cell boundaries in OpenGL coordinates 
            // The coordinate system ranges from -WIDTH to WIDTH, and -HEIGHT to HEIGHT.
            // A cell at (i, j) has a center at ((i * 2) - WIDTH, (j * 2) - HEIGHT)
            
            float x0 = i * 2.0f - GRID_W;
            float y0 = j * 2.0f - GRID_H;
            float x1 = (i + 1) * 2.0f - GRID_W;
            float y1 = (j + 1) * 2.0f - GRID_H;

            // Draw the quad for the cell
            glVertex2f(x0, y0);
            glVertex2f(x1, y0);
            glVertex2f(x1, y1);
            glVertex2f(x0, y1);
        }
    }
    glEnd();
}
void generateParticles(vector<vec3>& particles, int numParticles, int numClusters, float clusterRadius)
{
    particles.clear();

    // 1. Simulation bounds converted to 3D
    float minX = -engine.WIDTH;
    float maxX = engine.WIDTH;

    float minY = -engine.HEIGHT;
    float maxY = engine.HEIGHT;

    float minZ = -500;     // add this to your engine
    float maxZ = 500;

    // 2. Randomly pick cluster centers inside the 3D bounding box
    std::vector<vec3> clusterCenters;
    clusterCenters.reserve(numClusters);

    for (int i = 0; i < numClusters; ++i) {
        float x = minX + (std::rand() / (float)RAND_MAX) * (maxX - minX);
        float y = minY + (std::rand() / (float)RAND_MAX) * (maxY - minY);
        float z = minZ + (std::rand() / (float)RAND_MAX) * (maxZ - minZ);
        clusterCenters.emplace_back(x, y, z);
    }

    // 3. Global (uniform) vs clustered particles
    int globalParticles = numParticles / 3;
    int clusteredParticles = numParticles - globalParticles;

    // --- Global random 3D particles ---
    for (int i = 0; i < globalParticles; ++i) {
        float x = minX + (std::rand() / (float)RAND_MAX) * (maxX - minX);
        float y = minY + (std::rand() / (float)RAND_MAX) * (maxY - minY);
        float z = minZ + (std::rand() / (float)RAND_MAX) * (maxZ - minZ);
        particles.emplace_back(x, y, z);
    }

    // --- Clustered 3D particles (uniform spherical distribution) ---
    for (int i = 0; i < clusteredParticles; ++i) {
        const vec3& center = clusterCenters[std::rand() % numClusters];

        // Random direction on a sphere
        float u = (std::rand() / (float)RAND_MAX) * 2.0f - 1.0f; // cos(theta)
        float phi = (std::rand() / (float)RAND_MAX) * 2.0f * M_PI;

        float sqrt1MinusU2 = sqrt(1 - u * u);
        float dx = sqrt1MinusU2 * cos(phi);
        float dy = sqrt1MinusU2 * sin(phi);
        float dz = u;

        // Random radius inside sphere (cube-root preserves uniform density)
        float radius = clusterRadius * cbrt(std::rand() / (float)RAND_MAX);

        particles.emplace_back(
            center.x + dx * radius,
            center.y + dy * radius,
            center.z + dz * radius
        );
    }
}

void project_2d(const vector<vec3>& particles_3d, vector<vec2>& particles_2d) {
    particles_2d.clear();
    for (const vec3& p : particles_3d) {
        particles_2d.emplace_back(p.x, p.y); // simple orthographic projection onto XY plane
    }
}

// ================= Main ================= //
int main () {

    // ----- Generate particles -----
    vector<vec3> particles_3d;
    generateParticles(particles_3d, 10000, 5, 300.0f);

    vector<vec2> particles;
    project_2d(particles_3d, particles);
    
    // -- calculate density map --
    calculateDensityMap(particles);

    while (!glfwWindowShouldClose(engine.window)) {
        engine.run();

        // ---- 1. Draw the density map -------------------
        drawDensityMap();

        // ---- 2. draw particles on top ------------------
        glColor4f(1.0f, 1.0f, 1.0f, 0.3f);
        for (const vec2& p : particles)
            engine.drawFilledCircle(p.x, p.y, 1.0f);
        

        glfwSwapBuffers(engine.window);
        glfwPollEvents();
    }

    return 0;
}


