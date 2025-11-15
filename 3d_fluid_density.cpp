#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <thread>
#include <chrono>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
using namespace glm;
using namespace std;

const float c = 299792458.0f / 100000000.0f;    // speed of light in m/s
const float eu = 2.71828182845904523536f; // Euler's number
const float k = 8.9875517923e9f; // Coulomb's constant
const float a0 = 52.9f; // Bohr radius in pm
const float electron_r = 2.0f;
const float fieldRes = 25.0f;

// ================= Engine ================= //

struct Camera {
    // Center the camera orbit at (0, 0, 0) + init radius
    vec3 target = vec3(0.0f, 0.0f, 0.0f);
    float radius = 500.0f;

    // Camera angles
    float azimuth = 0.0f;
    float elevation = M_PI / 2.0f;

    // movement speeds
    float orbitSpeed = 0.01f;
    float panSpeed = 0.01f;
    double zoomSpeed = 125.0f;

    bool dragging = false;
    bool panning = false;
    bool moving = false; // For compute shader optimization
    double lastX = 0.0, lastY = 0.0;

    // Calculate camera position in world space
    vec3 position() const {
        float clampedElevation = clamp(elevation, 0.01f, float(M_PI) - 0.01f);
        // Orbit around (0,0,0) always
        return vec3(
            radius * sin(clampedElevation) * cos(azimuth),
            radius * cos(clampedElevation),
            radius * sin(clampedElevation) * sin(azimuth)
        );
    }
    void update() {
        // Always keep target at black hole center
        target = vec3(0.0f, 0.0f, 0.0f);
        if(dragging | panning) {
            moving = true;
        } else {
            moving = false;
        }
    }

    void processMouseMove(double x, double y) {
        float dx = float(x - lastX);
        float dy = float(y - lastY);

        if (dragging && panning) {
            // Pan: Shift + Left or Middle Mouse
            // Disable panning to keep camera centered on black hole
        }
        else if (dragging && !panning) {
            // Orbit: Left mouse only
            azimuth   += dx * orbitSpeed;
            elevation -= dy * orbitSpeed;
            elevation = clamp(elevation, 0.01f, float(M_PI) - 0.01f);
        }

        lastX = x;
        lastY = y;
        update();
    }
    void processMouseButton(int button, int action, int mods, GLFWwindow* win) {
        if (button == GLFW_MOUSE_BUTTON_LEFT || button == GLFW_MOUSE_BUTTON_MIDDLE) {
            if (action == GLFW_PRESS) {
                dragging = true;
                // Disable panning so camera always orbits center
                panning = false;
                glfwGetCursorPos(win, &lastX, &lastY);
            } else if (action == GLFW_RELEASE) {
                dragging = false;
                panning = false;
            }
        }
    }
    void processScroll(double xoffset, double yoffset) {
        radius -= yoffset * zoomSpeed;
        update();
    };
};
Camera camera;

struct Engine {

    // opengl vars
     GLFWwindow* window;
     GLuint shaderProgram;

    // vars - scale
    int WIDTH = 800;  // Window width
    int HEIGHT = 600; // Window height
    float width = 1000.0f; // Width of the viewport in picometers 
    float height = 750.0f; // Height of the viewport in picometers 
    

    Engine () {
        // init glfw
        if (!glfwInit()) {
            cerr << "GLFW init failed\n";
            exit(EXIT_FAILURE);
        }

        // create window
        window = glfwCreateWindow(WIDTH, HEIGHT, "Quantum Simulation by kavan G", nullptr, nullptr);
        if (!window) {
            cerr << "Failed to create GLFW window\n";
            glfwTerminate();
            exit(EXIT_FAILURE);
        }
        glViewport(0, 0, WIDTH, HEIGHT);
        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        glfwMakeContextCurrent(window);
        glEnable(GL_DEPTH_TEST);

        // Enable alpha blending for transparent objects
        glEnable(GL_BLEND);
        //glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE);

        // init glew
        glewExperimental = GL_TRUE;
        GLenum glewErr = glewInit();
        if (glewErr != GLEW_OK) {
            cerr << "Failed to initialize GLEW: "
                << (const char*)glewGetErrorString(glewErr)
                << "\n";
            glfwTerminate();
            exit(EXIT_FAILURE);
        }
        // Create shader program
        this->shaderProgram = CreateShaderProgram();
    }

    void run() {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glUseProgram(shaderProgram);

        mat4 view = lookAt(camera.position(), camera.target, vec3(0,1,0));
        mat4 projection = perspective(radians(45.0f), float(WIDTH)/HEIGHT, 0.1f, 10000.0f); // clipping distance

        glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "view"), 1, GL_FALSE, value_ptr(view));
        glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "projection"), 1, GL_FALSE, value_ptr(projection));

    };
    GLuint CreateShaderProgram(){
        const char* vertexShaderSource = R"(
        #version 330 core
        layout (location = 0) in vec3 aPos;

        uniform mat4 model;
        uniform mat4 view;
        uniform mat4 projection;
        uniform vec4 objectColor;

        out vec4 FragColorVS;

        void main() {
            gl_Position = projection * view * model * vec4(aPos, 1.0);
            FragColorVS = objectColor;
        })";

        const char* fragmentShaderSource = R"(
            #version 330 core
            in vec4 FragColorVS;
            out vec4 FragColor;

            void main() {
                FragColor = FragColorVS;
            }
        )";

        // vertex shader
        GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
        glShaderSource(vertexShader, 1, &vertexShaderSource, nullptr);
        glCompileShader(vertexShader);

        // fragment shader
        GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
        glShaderSource(fragmentShader, 1, &fragmentShaderSource, nullptr);
        glCompileShader(fragmentShader);

        GLuint shaderProgram = glCreateProgram();
        glAttachShader(shaderProgram, vertexShader);
        glAttachShader(shaderProgram, fragmentShader);
        glLinkProgram(shaderProgram);

        glDeleteShader(vertexShader);
        glDeleteShader(fragmentShader);

        return shaderProgram;
    };
    void CreateVBOVAO(GLuint& VAO, GLuint& VBO, const float* vertices, size_t vertexCount) {
        glGenVertexArrays(1, &VAO);
        glGenBuffers(1, &VBO);

        glBindVertexArray(VAO);
        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        glBufferData(GL_ARRAY_BUFFER, vertexCount * sizeof(float), vertices, GL_STATIC_DRAW);

        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);
        glBindVertexArray(0);
    }
    vec3 sphericalToCartesian(float r, float theta, float phi){
        float x = r * sin(theta) * cos(phi);
        float y = r * cos(theta);
        float z = r * sin(theta) * sin(phi);
        return vec3(x, y, z);
    };
    void SwitchTo2DRendering () {
        // --- SWITCH TO 2D ORTHO FOR DENSITY ---
        glm::mat4 ortho = glm::ortho(0.0f, float(WIDTH), 0.0f, float(HEIGHT));
        GLuint projLoc = glGetUniformLocation(shaderProgram, "projection");
        glUniformMatrix4fv(projLoc, 1, GL_FALSE, glm::value_ptr(ortho));

        glm::mat4 view2D = glm::mat4(1.0f);
        GLuint viewLoc = glGetUniformLocation(shaderProgram, "view");
        glUniformMatrix4fv(viewLoc, 1, GL_FALSE, glm::value_ptr(view2D));
    }
};
Engine engine;

void setupCameraCallbacks(GLFWwindow* window) {
    glfwSetWindowUserPointer(window, &camera);

    glfwSetMouseButtonCallback(window, [](GLFWwindow* win, int button, int action, int mods) {
        Camera* cam = (Camera*)glfwGetWindowUserPointer(win);
        cam->processMouseButton(button, action, mods, win);
    });

    glfwSetCursorPosCallback(window, [](GLFWwindow* win, double x, double y) {
        Camera* cam = (Camera*)glfwGetWindowUserPointer(win);
        cam->processMouseMove(x, y);
    });

    glfwSetScrollCallback(window, [](GLFWwindow* win, double xoffset, double yoffset) {
        Camera* cam = (Camera*)glfwGetWindowUserPointer(win);
        cam->processScroll(xoffset, yoffset);
    });
}


// ================= Objects ================= //
struct Grid {
    GLuint gridVAO, gridVBO;
    vector<float> vertices;
    Grid() {
        vertices = CreateGridVertices(500.0f, 2);
        engine.CreateVBOVAO(gridVAO, gridVBO, vertices.data(), vertices.size());
    }
    void Draw (GLint objectColorLoc) {
        glUseProgram(engine.shaderProgram);
        glUniform4f(objectColorLoc, 1.0f, 1.0f, 1.0f, 0.05f);
        glBindBuffer(GL_ARRAY_BUFFER, gridVBO);
        glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_DYNAMIC_DRAW);
        DrawGrid(engine.shaderProgram, gridVAO, vertices.size());
    }
    void DrawGrid(GLuint shaderProgram, GLuint gridVAO, size_t vertexCount) {
        glUseProgram(shaderProgram);
        glm::mat4 model = glm::mat4(1.0f); // Identity matrix for the grid
        GLint modelLoc = glGetUniformLocation(shaderProgram, "model");
        glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model));

        glBindVertexArray(gridVAO);
        glPointSize(2.0f);
        glDrawArrays(GL_LINES, 0, vertexCount / 3);
        glBindVertexArray(0);
    }
    vector<float> CreateGridVertices(float size, int divisions) {
        
        std::vector<float> vertices;
        float step = size / divisions;
        float halfSize = size / 2.0f;

        // amount to extend the central X-axis line (in same units as size)
        float extra = step * 3.0f; // adjust this factor to make the line stick out more/less
        int midZ = divisions / 2;

        // x axis
        for (int yStep = 3; yStep <= 3; ++yStep) {
            float y = 0;
            for (int zStep = 0; zStep <= divisions; ++zStep) {
                float z = -halfSize + zStep * step;
                for (int xStep = 0; xStep < divisions; ++xStep) {
                    float xStart = -halfSize + xStep * step;
                    float xEnd = xStart + step;

                    // If this is the central line (middle z), extend the very first and last segment
                    if (zStep == midZ) {
                        if (xStep == 0) {
                            xStart -= extra; // extend left end
                        }
                        if (xStep == divisions - 1) {
                            xEnd += extra;   // extend right end
                        }
                    }

                    vertices.push_back(xStart); vertices.push_back(y); vertices.push_back(z);
                    vertices.push_back(xEnd);   vertices.push_back(y); vertices.push_back(z);
                }
            }
        }
        // zaxis
        for (int xStep = 0; xStep <= divisions; ++xStep) {
            float x = -halfSize + xStep * step;
            for (int yStep = 3; yStep <= 3; ++yStep) {
                float y = 0;
                for (int zStep = 0; zStep < divisions; ++zStep) {
                    float zStart = -halfSize + zStep * step;
                    float zEnd = zStart + step;
                    vertices.push_back(x); vertices.push_back(y); vertices.push_back(zStart);
                    vertices.push_back(x); vertices.push_back(y); vertices.push_back(zEnd);
                }
            }
        }

        return vertices;

    }
};
Grid grid;

void drawParticle(vec3 particle, GLint modelLoc, GLint objectColorLoc) {
    // Draw each particle
    glUniform4f(objectColorLoc, 1.0f, 1.0f, 1.0f, 1.0f); // Red color for particles
    glm::mat4 model = glm::mat4(1.0f);
    model = glm::translate(model, particle);
    model = glm::scale(model, glm::vec3(2.0f)); // Scale to make it visible
    glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model));

    // Draw a simple point for the particle
    glPointSize(1.0f);
    glBegin(GL_POINTS);
    glVertex3f(0.0f, 0.0f, 0.0f);
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


const int GRID_W = 1000;
const int GRID_H = 1000;
const int GRID_D = 1000;
const float CELL_W = 100.0f;
const float CELL_H = 100.0f;
const float CELL_D = 100.0f;

vector<float> density(GRID_W * GRID_H * GRID_D, 0.0f);
float maxDensity = 1.0f;

// ------------------ Thermal colormap ------------------
vec3 densityToColour(float d, float maxD) {
    if (maxD <= 0.0f) maxD = 1.0f;
    float t = glm::clamp(d / maxD, 0.0f, 1.0f);
    float r, g, b;

    if (t < 0.33f) {
        float k = t / 0.33f;
        r = k; g = 0.0f; b = 0.0f;
    } else if (t < 0.66f) {
        float k = (t - 0.33f) / 0.33f;
        r = 1.0f; g = k; b = 0.0f;
    } else {
        float k = (t - 0.66f) / 0.34f;
        r = 1.0f; g = 1.0f; b = k;
    }

    return vec3(r, g, b);
}

// ------------------ Density Calculation ------------------
void calculateDensityMap(const vector<vec3>& particles) {
    std::fill(density.begin(), density.end(), 0.0f);
    maxDensity = 0.0f;

    const float INFLUENCE_RADIUS = 200.0f;
    const float INFLUENCE_RADIUS_SQ = INFLUENCE_RADIUS * INFLUENCE_RADIUS;
    const float DENSITY_SCALE = 1.0f;

    for (const vec3& p : particles) {
        int center_i = (int)round(p.x + GRID_W) / 2;
        int center_j = (int)round(p.y + GRID_H) / 2;
        int center_k = (int)round(p.z + GRID_D) / 2;

        int min_i = std::max(0, center_i - (int)ceil(INFLUENCE_RADIUS));
        int max_i = std::min(GRID_W - 1, center_i + (int)ceil(INFLUENCE_RADIUS));
        int min_j = std::max(0, center_j - (int)ceil(INFLUENCE_RADIUS));
        int max_j = std::min(GRID_H - 1, center_j + (int)ceil(INFLUENCE_RADIUS));
        int min_k = std::max(0, center_k - (int)ceil(INFLUENCE_RADIUS));
        int max_k = std::min(GRID_D - 1, center_k + (int)ceil(INFLUENCE_RADIUS));

        for (int k = min_k; k <= max_k; ++k) {
            for (int j = min_j; j <= max_j; ++j) {
                for (int i = min_i; i <= max_i; ++i) {
                    float cell_x = i * 2.0f - GRID_W;
                    float cell_y = j * 2.0f - GRID_H;
                    float cell_z = k * 2.0f - GRID_D;

                    float dx = p.x - cell_x;
                    float dy = p.y - cell_y;
                    float dz = p.z - cell_z;

                    float distSq = dx*dx + dy*dy + dz*dz;

                    if (distSq < INFLUENCE_RADIUS_SQ) {
                        float dist = sqrt(distSq);
                        float influence = (INFLUENCE_RADIUS - dist) / INFLUENCE_RADIUS;
                        float kernel_weight = influence * influence * DENSITY_SCALE;

                        int index = k * GRID_W * GRID_H + j * GRID_W + i;
                        density[index] += kernel_weight;
                        maxDensity = std::max(maxDensity, density[index]);
                    }
                }
            }
        }
    }
}

// ------------------ Rendering 3D Density Map ------------------
void drawDensityMap() {
    glBegin(GL_QUADS); // For simplicity, you can later switch to 3D cubes (GL_QUADS + z)
    for (int k = 0; k < GRID_D; ++k) {
        for (int j = 0; j < GRID_H; ++j) {
            for (int i = 0; i < GRID_W; ++i) {
                int index = k * GRID_W * GRID_H + j * GRID_W + i;
                float d = density[index];
                vec3 color = densityToColour(d, maxDensity);
                glColor3f(color.r, color.g, color.b);

                // You can draw a small cube/quad at each grid cell
                float x0 = i * 2.0f - GRID_W;
                float y0 = j * 2.0f - GRID_H;
                float z0 = k * 2.0f - GRID_D;
                float x1 = x0 + CELL_W;
                float y1 = y0 + CELL_H;
                float z1 = z0 + CELL_D;

                // Draw one face (XY) for now; extend for full cube if needed
                glVertex3f(x0, y0, z0);
                glVertex3f(x1, y0, z0);
                glVertex3f(x1, y1, z0);
                glVertex3f(x0, y1, z0);
            }
        }
    }
    glEnd();
}

void drawElectronCloud(const vector<vec3>& electrons, GLint modelLoc, GLint objectColorLoc) {
    for (const vec3& e : electrons) {
        float alpha = 1.0f; // make them transparent
        glUniform4f(objectColorLoc, 0.0f, 0.5f, 1.0f, alpha);

        glm::mat4 model = glm::mat4(1.0f);
        model = glm::translate(model, e);
        model = glm::scale(model, glm::vec3(20.0f));
        glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model));

        glBegin(GL_POINTS);
        glVertex3f(0.0f, 0.0f, 0.0f);
        glEnd();
    }
}

void project_particles(const vector<vec3>& particles_3d, vector<vec2>& particles_2d) {
    particles_2d.clear();
    particles_2d.reserve(particles_3d.size());

    for (const vec3& p : particles_3d) {
        particles_2d.emplace_back(p.x, p.y); // Simple orthographic projection
    }
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
float radialProbability3p(float r) {
    float x = r / a0;
    float R = x * (1.0f - x / 6.0f) * exp(-x / 3.0f);
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

void sample1s(vector<vec3> &particles) {   // change return type to void
    float r = sampleR1s();
    float theta = acos(1.0f - 2.0f * static_cast<float>(rand()) / RAND_MAX); // [0, pi]
    float phi   = 2.0f * M_PI * static_cast<float>(rand()) / RAND_MAX;       // [0, 2pi]
    vec3 electronPos = engine.sphericalToCartesian(r, theta, phi);
    
    // Construct Particle in-place
    if (true) {
        particles.emplace_back(  electronPos ); // cyan-ish color
    }
}
void sample2s(vector<vec3> &particles) {
    float r = sampleR2s();
    float theta = acos(1.0f - 2.0f * static_cast<float>(rand()) / RAND_MAX); // [0, π]
    float phi   = 2.0f * M_PI * static_cast<float>(rand()) / RAND_MAX;       // [0, 2π]
    
    vec3 electronPos = engine.sphericalToCartesian(r, theta, phi);

    // Use a distinct color for visualization, e.g. yellow for 2s
    if (true) {
    particles.emplace_back( electronPos); // cyan-ish color
    }
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

    vec3 electronPos = engine.sphericalToCartesian(r, theta, phi);

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

    vec3 electronPos = engine.sphericalToCartesian(r, theta, phi);

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

    vec3 electronPos = engine.sphericalToCartesian(r, theta, phi);

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

    vec3 electronPos = engine.sphericalToCartesian(r, theta, phi);

    // visualize (different color for 3p_z)
    particles.emplace_back( electronPos ); // cyan-ish for 3p_z
}

// ================= Main ================= //
int main () {
    setupCameraCallbacks(engine.window);
    GLint modelLoc = glGetUniformLocation(engine.shaderProgram, "model");
    GLint objectColorLoc = glGetUniformLocation(engine.shaderProgram, "objectColor");
    glUseProgram(engine.shaderProgram);

    vector<vec3> particles;
                // list dumbass  num    clus   rad
    //generateParticles(particles, 5000, 5, 500.0f);

    vector<vec2> particles_projection;


    //calculateDensityMap(particles);
    for (int i = 0; i < 2000; ++i) {
        //sample1s(particles);
        //sample2s();
        //sample2p_x();
        sample2p_y(particles);
        //sample2p_z();
        sample3p_z(particles);
        sample3p_z(particles);
        sample3p_z(particles);
        sample3p_z(particles);
        sample3p_z(particles);
        sample3p_z(particles);
    }

    // ------------------ RENDERING LOOP ------------------
    while (!glfwWindowShouldClose(engine.window)) {
        engine.run();

        // ---- DRAW GRID ----
        grid.Draw(objectColorLoc);

        // ---- 1. Draw the density map -------------------
        //drawDensityMap();
        drawElectronCloud(particles, modelLoc, objectColorLoc);


        // ---- 2. Draw particles -------------------
        for (const auto& particle : particles) {
            drawParticle(particle, modelLoc, objectColorLoc);
        }

        glfwSwapBuffers(engine.window);
        glfwPollEvents();
    }
    
    glfwDestroyWindow(engine.window);
    glfwTerminate();
    return 0;
}