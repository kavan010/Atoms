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

const float a0 = 52.9f; // Bohr radius in pm

struct Camera {
    vec3 target = vec3(0.0f);
    vec3 up = vec3(0,1,0);

    // make radius positive and on a sensible scale so zoom is meaningful
    float radius = 200.0f;

    float azimuth = 0.0f;
    float elevation = M_PI / 2.0f;

    float orbitSpeed = 0.01f;
    double zoomFactor = 2.0;

    bool dragging = false;
    double lastX = 0.0, lastY = 0.0;

    vec3 getPosition() const {
        float x = radius * sin(elevation) * cos(azimuth);
        float y = radius * cos(elevation);
        float z = radius * sin(elevation) * sin(azimuth);
        return target + vec3(x, y, z);
    }

    void rotate(float dAzimuth, float dElevation) {
        azimuth += dAzimuth;
        elevation += dElevation;

        if(elevation < 0.01f) elevation = 0.01f;
        if(elevation > M_PI - 0.01f) elevation = M_PI - 0.01f;

        if(azimuth > 2*M_PI) azimuth -= 2*M_PI;
        if(azimuth < 0) azimuth += 2*M_PI;
    }

    void processMouseMove(double x, double y) {
        if(dragging) {
            double dx = x - lastX;
            double dy = y - lastY;
            rotate(-dx * orbitSpeed, -dy * orbitSpeed); // fixed axes
        }
        lastX = x;
        lastY = y;
    }
    void processMouseButton(int button, int action, int mods, GLFWwindow* win) {
        if(button == GLFW_MOUSE_BUTTON_LEFT) {
            if(action == GLFW_PRESS) {
                dragging = true;
                // get cursor in GLFW coords (origin top-left) then convert to OpenGL coords (origin bottom-left)
                double lx, ly;
                glfwGetCursorPos(win, &lx, &ly);
                int ww, wh;
                glfwGetWindowSize(win, &ww, &wh);
                lastX = lx;
                lastY = (double)wh - ly;
            } else if(action == GLFW_RELEASE) {
                dragging = false;
            }
        }
    }

    void processScroll(double xoffset, double yoffset) {
       // multiplicative zoom — smaller step for finer control
       const float zoomStep = 0.10f; // 10% per notch
       float factor;
       if (yoffset > 0.0) factor = powf(1.0f - zoomStep, (float) yoffset);   // zoom in
       else factor = powf(1.0f + zoomStep, (float) -yoffset);               // zoom out
       radius *= factor;

       if (radius < 1.0f) radius = 1.0f;
       if (radius > 50000.0f) radius = 50000.0f;
   }
    vec2 projectTo2D(const vec3& p, float left, float right, float bottom, float top) const {
        vec3 camPos = getPosition();
        mat4 view = lookAt(camPos, target, up);
        mat4 proj = ortho(left, right, bottom, top, -1.0f, 1.0f);

        vec4 clip = proj * view * vec4(p, 1.0f);

        // NDC
        float xN = clip.x / clip.w;
        float yN = clip.y / clip.w;

        // flip Y so result uses bottom-left origin like your drawing code
        yN = -yN;

        // Scale to your engine units (match engine.run projection below)
        float x = left + (xN + 1.0f) * 0.5f * (right - left);
        float y = bottom + (yN + 1.0f) * 0.5f * (top - bottom);

        return vec2(x, y);
    }
};
Camera camera;

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
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // Use camera.radius to control orthographic extents (zoom)
        // Choose a base scale so radius==200 gives the original extents

        const float baseRadius = 200.0f;
        float zoomScale = camera.radius / baseRadius;
        if (zoomScale <= 0.001f) zoomScale = 0.001f;
        double left   = -WIDTH * zoomScale;
        double right  =  WIDTH * zoomScale;
        double bottom = -HEIGHT * zoomScale;
        double top    =  HEIGHT * zoomScale;

        // Save current matrices
        glMatrixMode(GL_PROJECTION);
        glPushMatrix();
        glLoadIdentity();
        glOrtho(left, right, bottom, top, -1.0, 1.0);
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

    vec3 sphericalToCartesian(float r, float theta, float phi){
        float x = r * sin(theta) * cos(phi);
        float y = r * cos(theta);
        float z = r * sin(theta) * sin(phi);
        return vec3(x, y, z);
    };
};
Engine engine;


struct Atom {

    vec3 pos;
    Atom(vec3 pos) {
        this->pos = pos;
    }
    // ================= Probability functions ================= //
    float radialProbability1s(float r) {
        float r_bohr = r / a0;
        return 4.0 * r_bohr*r_bohr * exp(-2.0 * r_bohr);
        //4 * r**2 * np.exp(-2 * r)
    }
    float radialProbability2s(float r) {
        float r_bohr = r / a0;  // convert to Bohr radii
        return 0.5f * pow(r_bohr, 2) * pow(1 - r_bohr / 2.0f, 2) * exp(-r_bohr);
    }
    float radialProbability2p(float r) {
        float r_bohr = r / a0;  // convert to Bohr radii
        return (pow(r_bohr, 4) / 24.0f) * exp(-r_bohr);
    }
    float radialProbability3s(float r) {
        float x = r / a0;
        float R = (1.0f - (2.0f/3.0f)*x + (2.0f/27.0f)*x*x) * exp(-x / 3.0f);
        return R * R * r * r;
    }
    float radialProbability3p(float r) {
        float x = r / a0;
        float R = x * (1.0f - x / 6.0f) * exp(-x / 3.0f);
        return R * R * r * r; // multiply by r^2 for probability density
    }
    float radialProbability3dxy(float r) {
        float x = r / a0;
        float R = x*x * (6.0f - x) * exp(-x / 3.0f); // 3d radial part
        return R * R * r * r; // multiply by r^2 for probability density
    }

    float sampleR1s() {
        float r_max = 5.0f * a0;  // arbitrary max radius
        float P_max = radialProbability1s(0.0f); // max occurs at r ~ a0 (approx)
        
        while (true) {
            float r = static_cast<float>(rand()) / RAND_MAX * r_max;
            float y = static_cast<float>(rand()) / RAND_MAX * P_max;
            if (y <= radialProbability1s(r)) return r;
        }
    }
    float sampleR2s() {
        float r_max = 10.0f * a0;  // 2s extends farther out than 1s
        float P_max = radialProbability2s(a0); // rough peak near r ≈ a0

        while (true) {
            float r = static_cast<float>(rand()) / RAND_MAX * r_max;
            float y = static_cast<float>(rand()) / RAND_MAX * P_max;
            if (y <= radialProbability2s(r)) return r;
        }
    }
    float sampleR3s() {
        float r_max = 20.0f * a0;  // 3s extends farther than 1s/2s
        float P_max = radialProbability3s(3.0f * a0); // near the 3s radial peak

        while (true) {
            float r = static_cast<float>(rand()) / RAND_MAX * r_max;
            float y = static_cast<float>(rand()) / RAND_MAX * P_max;
            if (y <= radialProbability3s(r)) return r;
        }
    }
    float sampleR2p() {
        float r_max = 15.0f * a0;
        float P_max = radialProbability2p(4.0f * a0);

        while (true) {
            float r = static_cast<float>(rand()) / RAND_MAX * r_max;
            float y = static_cast<float>(rand()) / RAND_MAX * P_max;
            if (y <= radialProbability2p(r)) return r;
        }
    }
    float sampleR3p() {
        float r_max = 25.0f * a0; // 3p orbitals extend farther than 2p
        float P_max = radialProbability3p(8.0f * a0); // estimate peak around ~8a0

        while (true) {
            float r = static_cast<float>(rand()) / RAND_MAX * r_max;
            float y = static_cast<float>(rand()) / RAND_MAX * P_max;
            if (y <= radialProbability3p(r)) return r;
        }
    }
    float sampleR3dxy() {
        float r_max = 30.0f * a0; // slightly bigger than 3p
        float P_max = radialProbability3dxy(12.0f * a0); // rough peak estimate

        while (true) {
            float r = static_cast<float>(rand()) / RAND_MAX * r_max;
            float y = static_cast<float>(rand()) / RAND_MAX * P_max;
            if (y <= radialProbability3dxy(r)) return r;
        }
    }

    void sample1s(vector<vec3> &particles) {   // change return type to void
        float r = sampleR1s();
        float theta = acos(1.0f - 2.0f * static_cast<float>(rand()) / RAND_MAX); // [0, pi]
        float phi   = 2.0f * M_PI * static_cast<float>(rand()) / RAND_MAX;       // [0, 2pi]
        vec3 electronPos = engine.sphericalToCartesian(r, theta, phi) + pos;
        
        // Construct Particle in-place
        if (true) {
            particles.emplace_back(  electronPos ); // cyan-ish color
        }
    }
    void sample2s(vector<vec3> &particles) {
        float r = sampleR2s();
        float theta = acos(1.0f - 2.0f * static_cast<float>(rand()) / RAND_MAX); // [0, π]
        float phi   = 2.0f * M_PI * static_cast<float>(rand()) / RAND_MAX;       // [0, 2π]
        
        vec3 electronPos = engine.sphericalToCartesian(r, theta, phi) + pos;

        // Use a distinct color for visualization, e.g. yellow for 2s
        if (true) {
        particles.emplace_back( electronPos); // cyan-ish color
        }
    }
    void sample3s(vector<vec3> &particles) {
        float r = sampleR3s();

        // Uniform spherical angle distribution
        float theta = acos(1.0f - 2.0f * static_cast<float>(rand()) / RAND_MAX);
        float phi   = 2.0f * M_PI * static_cast<float>(rand()) / RAND_MAX;

        vec3 pos = engine.sphericalToCartesian(r, theta, phi) + pos;

        particles.emplace_back(pos);
    }
    void sample2p_x(vector<vec3> &particles) {
        float r = sampleR2p();

        float theta, phi;

        // --- sample theta ---
        while (true) {
            theta = acos(1.0f - 2.0f * static_cast<float>(rand()) / RAND_MAX); // [0, pi]
            float prob = pow(sin(theta), 3); // sin^3(theta)
            if (static_cast<float>(rand()) / RAND_MAX <= prob) break;
        }

        // --- sample phi ---
        while (true) {
            phi = 2.0f * M_PI * static_cast<float>(rand()) / RAND_MAX; // [0, 2pi]
            float prob = pow(cos(phi), 2); // cos^2(phi)
            if (static_cast<float>(rand()) / RAND_MAX <= prob) break;
        }

        vec3 electronPos = engine.sphericalToCartesian(r, theta, phi) + pos;

        // Keep one lobe for visualization
        if (true) {
            particles.emplace_back(  electronPos ); // red for 2p_x
        }
    }
    void sample2p_y(vector<vec3> &particles) {
        float r = sampleR2p();

        float theta, phi;

        // --- sample theta ---
        while (true) {
            theta = acos(1.0f - 2.0f * static_cast<float>(rand()) / RAND_MAX); // [0, pi]
            float prob = pow(sin(theta), 3); // sin^3(theta)
            if (static_cast<float>(rand()) / RAND_MAX <= prob) break;
        }

        // --- sample phi --- (changed to sin^2 to orient along Y axis)
        while (true) {
            phi = 2.0f * M_PI * static_cast<float>(rand()) / RAND_MAX; // [0, 2pi]
            float prob = pow(sin(phi), 2); // sin^2(phi) -> aligns lobes along Y
            if (static_cast<float>(rand()) / RAND_MAX <= prob) break;
        }

        vec3 electronPos = engine.sphericalToCartesian(r, theta, phi) + pos;

        // Keep one lobe for visualization
        if (true) {
            particles.emplace_back( electronPos ); // red for 2p_y (keeps existing color)
        }
    }
    void sample2p_z(vector<vec3> &particles) {
        float r = sampleR2p();
        float theta, phi;

        // --- sample theta --- (weight with cos^2 to align lobes along Z)
        while (true) {
            theta = acos(1.0f - 2.0f * static_cast<float>(rand()) / RAND_MAX); // [0, pi]
            float prob = pow(cos(theta), 2); // cos^2(theta)
            if (static_cast<float>(rand()) / RAND_MAX <= prob) break;
        }

        // --- sample phi --- (uniform)
        phi = 2.0f * M_PI * static_cast<float>(rand()) / RAND_MAX; // [0, 2pi]

        vec3 electronPos = engine.sphericalToCartesian(r, theta, phi) + pos;

        // Keep one lobe for visualization
        particles.emplace_back(electronPos); // red for 2p_z
    }

    void sample3p_z(vector<vec3> &particles) {
        float r = sampleR3p(); // 3p radial distribution
        float theta, phi;

        // --- sample theta --- (same angular part as 2p_z)
        while (true) {
            theta = acos(1.0f - 2.0f * static_cast<float>(rand()) / RAND_MAX); // [0, pi]
            float prob = pow(cos(theta), 2); // cos^2(theta) for p_z alignment
            if (static_cast<float>(rand()) / RAND_MAX <= prob) break;
        }

        // --- sample phi --- (uniform)
        phi = 2.0f * M_PI * static_cast<float>(rand()) / RAND_MAX; // [0, 2pi]

        vec3 electronPos = engine.sphericalToCartesian(r, theta, phi) + pos;

        // visualize (different color for 3p_z)
        particles.emplace_back( electronPos ); // cyan-ish for 3p_z
    }
    void sample3dxy(vector<vec3> &particles) {
        float r = sampleR3dxy();
        float theta, phi;

        while (true) {
            theta = acos(1.0f - 2.0f * static_cast<float>(rand()) / RAND_MAX); // [0, pi]
            phi = 2.0f * M_PI * static_cast<float>(rand()) / RAND_MAX;         // [0, 2pi]

            float prob = pow(sin(theta), 2) * pow(sin(phi), 2); // 3dxy angular distribution
            if (static_cast<float>(rand()) / RAND_MAX <= prob) break;
        }

        vec3 electronPos = engine.sphericalToCartesian(r, theta, phi) + pos;
        particles.emplace_back(electronPos); // can color differently for 3dxy
    }
};
vector<Atom> atoms = 
{
    Atom(vec3(0.0f, 0.0f, 0.0f)),
    Atom(vec3(0.0f, 2500.0f, 0.0f))
};

// ================= Callbacks ================= //
void setupCameraCallbacks(GLFWwindow* window) {
    glfwSetWindowUserPointer(window, &camera);

    glfwSetMouseButtonCallback(window, [](GLFWwindow* win, int button, int action, int mods) {
        Camera* cam = (Camera*)glfwGetWindowUserPointer(win);
        cam->processMouseButton(button, action, mods, win);
    });

    // flip Y here so Camera sees bottom-left origin coordinates
    glfwSetCursorPosCallback(window, [](GLFWwindow* win, double x, double y) {
        Camera* cam = (Camera*)glfwGetWindowUserPointer(win);
        int ww, wh;
        glfwGetWindowSize(win, &ww, &wh);
        double yf = (double)wh - y;
        cam->processMouseMove(x, yf);
    });

    glfwSetScrollCallback(window, [](GLFWwindow* win, double xoffset, double yoffset) {
        Camera* cam = (Camera*)glfwGetWindowUserPointer(win);
        cam->processScroll(xoffset, yoffset);
    });
}


// ================= Grid Vars ================= //
const int   GRID_W = engine.WIDTH;
const int   GRID_H = engine.HEIGHT;
const float CELL_W = 1.0f;
const float CELL_H = 1.0f;

vector<float> density(GRID_W * GRID_H, 0.0f);
float maxDensity = 1.0f;

// ================= Density map functions ================= //
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
    // Reset
    fill(density.begin(), density.end(), 0.0f);
    maxDensity = 0.0f;

    // Match ortho extents used by engine.run
    const float baseRadius = 200.0f;
    float zoomScale = camera.radius / baseRadius;
    if (zoomScale <= 0.001f) zoomScale = 0.001f;
    float left   = -engine.WIDTH * zoomScale;
    float right  =  engine.WIDTH * zoomScale;
    float bottom = -engine.HEIGHT * zoomScale;
    float top    =  engine.HEIGHT * zoomScale;

    // Cell size in world coordinates
    float cellW = (right - left) / float(GRID_W);
    float cellH = (top - bottom) / float(GRID_H);

    // Influence radius in world units (tweak if you want wider/narrower blur)
    const float INFLUENCE_RADIUS = 200.0f;
    const float INFLUENCE_RADIUS_SQ = INFLUENCE_RADIUS * INFLUENCE_RADIUS;

    for (const vec2& p : particles) {
        // particle outside current extents -> skip
        if (p.x < left || p.x > right || p.y < bottom || p.y > top) continue;

        int center_i = int(floor((p.x - left) / cellW));
        int center_j = int(floor((p.y - bottom) / cellH));
        if (center_i < 0 || center_i >= GRID_W || center_j < 0 || center_j >= GRID_H) continue;

        int radiusCellsX = (int)ceil(INFLUENCE_RADIUS / cellW);
        int radiusCellsY = (int)ceil(INFLUENCE_RADIUS / cellH);

        int min_i = std::max(0, center_i - radiusCellsX);
        int max_i = std::min(GRID_W - 1, center_i + radiusCellsX);
        int min_j = std::max(0, center_j - radiusCellsY);
        int max_j = std::min(GRID_H - 1, center_j + radiusCellsY);

        for (int j = min_j; j <= max_j; ++j) {
            for (int i = min_i; i <= max_i; ++i) {
                // cell center in world coords
                float cell_x = left + (i + 0.5f) * cellW;
                float cell_y = bottom + (j + 0.5f) * cellH;

                float dx = p.x - cell_x;
                float dy = p.y - cell_y;
                float distSq = dx*dx + dy*dy;

                if (distSq < INFLUENCE_RADIUS_SQ) {
                    float dist = sqrtf(distSq);
                    float influence = (INFLUENCE_RADIUS - dist) / INFLUENCE_RADIUS;
                    float kernel_weight = influence * influence;
                    int index = j * GRID_W + i;
                    density[index] += kernel_weight;
                    maxDensity = std::max(maxDensity, density[index]);
                }
            }
        }
    }
}
void drawDensityMap() {
    // Match ortho extents used by engine.run
    const float baseRadius = 200.0f;
    float zoomScale = camera.radius / baseRadius;
    if (zoomScale <= 0.001f) zoomScale = 0.001f;
    float left   = -engine.WIDTH * zoomScale;
    float right  =  engine.WIDTH * zoomScale;
    float bottom = -engine.HEIGHT * zoomScale;
    float top    =  engine.HEIGHT * zoomScale;

    float cellW = (right - left) / float(GRID_W);
    float cellH = (top - bottom) / float(GRID_H);

    glBegin(GL_QUADS);
    for (int j = 0; j < GRID_H; ++j) {
        for (int i = 0; i < GRID_W; ++i) {
            int index = j * GRID_W + i;
            float d = density[index];
            vec3 color = densityToColour(d, maxDensity);
            glColor3f(color.r, color.g, color.b);

            float x0 = left + i * cellW;
            float y0 = bottom + j * cellH;
            float x1 = x0 + cellW;
            float y1 = y0 + cellH;

            glVertex2f(x0, y0);
            glVertex2f(x1, y0);
            glVertex2f(x1, y1);
            glVertex2f(x0, y1);
        }
    }
    glEnd();
}

void project_2d(const vector<vec3>& particles_3d, vector<vec2>& particles_2d) {
    particles_2d.clear();

    // match the ortho extents used in Engine::run so projection/viewport align
    const float baseRadius = 200.0f;
    float zoomScale = camera.radius / baseRadius;
    if (zoomScale <= 0.001f) zoomScale = 0.001f;

    float left   = -engine.WIDTH * zoomScale;
    float right  =  engine.WIDTH * zoomScale;
    float bottom = -engine.HEIGHT * zoomScale;
    float top    =  engine.HEIGHT * zoomScale;

    for (const vec3& p : particles_3d) {
        particles_2d.push_back(camera.projectTo2D(p, left, right, bottom, top));
    }
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

// ================= Main ================= //
int main () {
    setupCameraCallbacks(engine.window);



    // ----- Generate particles -----
    vector<vec3> particles_3d;
    //generateParticles(particles_3d, 10000, 5, 300.0f);
    for (Atom & atom : atoms) {
        for (int i = 0; i < 2000; i++) {
            //atom.sample1s(particles_3d);
            // atom.sample2p_y(particles_3d);
            // atom.sample3s(particles_3d);
            atom.sample3p_z(particles_3d);
            // atom.sample3dxy(particles_3d);
        }
    }

    // ------- 1. Declare 2D projection -------------------
    vector<vec2> particles;

    while (!glfwWindowShouldClose(engine.window)) {
        engine.run();

        // particles_3d.clear();
        // for (Atom & atom : atoms) {
        //     atom.sample3p_z(particles_3d);
        // }

        // -------- project 3D to 2D and recalculate density --------
        project_2d(particles_3d, particles);
        calculateDensityMap(particles);

        

        // ------- 1. Draw the density map -------------------
        drawDensityMap();

        // ------- 2. Draw the particles -------------------
 
        // Draw 2D projected particles in this orthographic space
        // glColor4f(1.0f, 1.0f, 1.0f, 0.7f);
        // for (const vec2& p : particles)
        //     engine.drawFilledCircle(p.x, p.y, 1.0f);

        glfwSwapBuffers(engine.window);
        glfwPollEvents();
    }

    return 0;
}


