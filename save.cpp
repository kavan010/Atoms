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
#include <iomanip>
#include <thread>
#include <chrono>
#include "../libs/QuickHull/QuickHull.hpp"
#include <unordered_map>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
using namespace glm;
using namespace std;

const float c = 299792458.0f / 100000000.0f;    // speed of light in m/s
const float eu = 2.71828182845904523536f; // Euler's number
const float k = 8.9875517923e9f; // Coulomb's constant
const float a0 = 52.9f; // Bohr radius in pm
const float electron_r = 1.0f;
const float fieldRes = 10.0f;

// ================= Engine ================= //
struct Particle;

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
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

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
        mat4 projection = perspective(radians(45.0f), float(WIDTH)/HEIGHT, 0.1f, 2000.0f);

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

        // x axis
        for (int yStep = 3; yStep <= 3; ++yStep) {
            float y = 0;
            for (int zStep = 0; zStep <= divisions; ++zStep) {
                float z = -halfSize + zStep * step;
                for (int xStep = 0; xStep < divisions; ++xStep) {
                    float xStart = -halfSize + xStep * step;
                    float xEnd = xStart + step;
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

struct Particle {
    GLuint VAO, VBO;
    float radius; // in femtometers
    vec4 color; // RGB values between 0 and 1
    vec3 position; // in picometers
    vector<float> vertices = DrawSphere();

    Particle(float r, vec4 col, vec3 pos) : radius(r), color(col), position(pos) {
        engine.CreateVBOVAO(VAO, VBO, vertices.data(), vertices.size());
    }
    void Draw (GLint objectColorLoc, GLint modelLoc) {
        glUniform4f(objectColorLoc, color.r, color.g, color.b, color.a);
        glm::mat4 model = glm::mat4(1.0f);
        model = glm::translate(model, position); // Apply position here
        glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model));
        glUniform1i(glGetUniformLocation(engine.shaderProgram, "GLOW"), 0);
        glBindVertexArray(VAO);
    }
    vector<float> DrawSphere() {
        std::vector<float> vertices;
        int stacks = 25;
        int sectors = 25;

        for(float i = 0.0f; i <= stacks; ++i){
            float theta1 = (i / stacks) * pi<float>();
            float theta2 = (i+1) / stacks * pi<float>();
            for (float j = 0.0f; j < sectors; ++j){
                float phi1 = j / sectors * 2 * glm::pi<float>();
                float phi2 = (j+1) / sectors * 2 * glm::pi<float>();
                glm::vec3 v1 = engine.sphericalToCartesian(this->radius, theta1, phi1);
                glm::vec3 v2 = engine.sphericalToCartesian(this->radius, theta1, phi2);
                glm::vec3 v3 = engine.sphericalToCartesian(this->radius, theta2, phi1);
                glm::vec3 v4 = engine.sphericalToCartesian(this->radius, theta2, phi2);    
                // Triangle 1: v1-v2-v3
                vertices.insert(vertices.end(), {v1.x, v1.y, v1.z}); //      /|
                vertices.insert(vertices.end(), {v2.x, v2.y, v2.z}); //     / |
                vertices.insert(vertices.end(), {v3.x, v3.y, v3.z}); //    /__|

                // Triangle 2: v2-v4-v3
                vertices.insert(vertices.end(), {v2.x, v2.y, v2.z});
                vertices.insert(vertices.end(), {v4.x, v4.y, v4.z});
                vertices.insert(vertices.end(), {v3.x, v3.y, v3.z});
            }   
        }
        return vertices;
    }
};
vector<Particle> particles{
            // r   // color                      // position
    Particle(8.7f, vec4(1.0f, 0.0f, 0.0f, 1.0f), vec3(0.0f, 0.0f, 0.0f)), // nucleus
};

// --- ADD: simple Cube type for red cubes ---
struct Cube {
    GLuint VAO = 0u, VBO = 0u;
    vec3 position;
    vec4 color;
    vector<float> vertices;
    float size = 20.0f;        // store size so we can compute bounds
    int particleCount = 0;     // current cached count

    Cube(vec3 pos, vec4 col, float size_ = 20.0f) : position(pos), color(col), size(size_) {
        float hs = size * 0.5f;
        // 12 triangles (36 vertices) for a cube centered at origin
        vertices = {
            -hs,-hs,-hs,  hs,-hs,-hs,  hs, hs,-hs,
             hs, hs,-hs, -hs, hs,-hs, -hs,-hs,-hs,

            -hs,-hs, hs,  hs,-hs, hs,  hs, hs, hs,
             hs, hs, hs, -hs, hs, hs, -hs,-hs, hs,

            -hs, hs,-hs, -hs, hs, hs, -hs,-hs, hs,
            -hs,-hs, hs, -hs,-hs,-hs, -hs, hs,-hs,

             hs, hs,-hs,  hs, hs, hs,  hs,-hs, hs,
             hs,-hs, hs,  hs,-hs,-hs,  hs, hs,-hs,

            -hs,-hs,-hs, -hs,-hs, hs,  hs,-hs, hs,
             hs,-hs, hs,  hs,-hs,-hs, -hs,-hs,-hs,

            -hs, hs,-hs,  hs, hs,-hs,  hs, hs, hs,
             hs, hs, hs, -hs, hs, hs, -hs, hs,-hs
        };
        engine.CreateVBOVAO(VAO, VBO, vertices.data(), vertices.size());
    }

    void Draw(GLint objectColorLoc, GLint modelLoc) {
        glUniform4f(objectColorLoc, color.r, color.g, color.b, color.a);
        glm::mat4 model = glm::mat4(1.0f);
        model = glm::translate(model, position);
        glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model));
        glBindVertexArray(VAO);
        glDrawArrays(GL_TRIANGLES, 0, (GLsizei)(vertices.size() / 3));
        glBindVertexArray(0);
    }
};

// ADD global cubes container
vector<Cube> cubes;


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

void sample1s() {   // change return type to void
    float r = sampleR1s();
    float theta = acos(1.0f - 2.0f * static_cast<float>(rand()) / RAND_MAX); // [0, pi]
    float phi   = 2.0f * M_PI * static_cast<float>(rand()) / RAND_MAX;       // [0, 2pi]
    vec3 electronPos = engine.sphericalToCartesian(r, theta, phi);
    
    // Construct Particle in-place
    if (electronPos.z > 0 || electronPos.y < 0) {
    particles.emplace_back(electron_r, vec4(0.0f, 1.0f, 1.0f, 1.0f), electronPos); // cyan-ish color
    }
}
void sample2s() {
    float r = sampleR2s();
    float theta = acos(1.0f - 2.0f * static_cast<float>(rand()) / RAND_MAX); // [0, π]
    float phi   = 2.0f * M_PI * static_cast<float>(rand()) / RAND_MAX;       // [0, 2π]
    
    vec3 electronPos = engine.sphericalToCartesian(r, theta, phi);

    // Use a distinct color for visualization, e.g. yellow for 2s
    if (electronPos.z > 0 || electronPos.y < 0) {
    particles.emplace_back(electron_r, vec4(0.0f, 1.0f, 1.0f, 1.0f), electronPos); // cyan-ish color
    }
}
void sample2p_x() {
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
    if (electronPos.z > 0 || electronPos.y < 0) {
        particles.emplace_back(electron_r, vec4(1.0f, 0.0f, 0.0f, 1.0f), electronPos); // red for 2p_x
    }
}
void sample2p_y() {
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
        particles.emplace_back(electron_r, vec4(1.0f, 0.0f, 1.0f, 1.0f), electronPos); // red for 2p_y (keeps existing color)
    }
}
void sample2p_z() {
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
    particles.emplace_back(electron_r, vec4(0.0f, 1.0f, 1.0f, 1.0f), electronPos); // red for 2p_z
}

void sample3p_z() {
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
    particles.emplace_back(electron_r, vec4(0.0f, 1.0f, 1.0f, 1.0f), electronPos); // cyan-ish for 3p_z
}
// ================= Main ================= //
int main () {
    setupCameraCallbacks(engine.window);
    GLint modelLoc = glGetUniformLocation(engine.shaderProgram, "model");
    GLint objectColorLoc = glGetUniformLocation(engine.shaderProgram, "objectColor");
    glUseProgram(engine.shaderProgram);

    // ---- GENERATE PARTICLES ---- //
    for (int i = 0; i < 2000; ++i) {
        //sample1s();
        //sample2s();
        //sample2p_x();
        //sample2p_y();
        //sample2p_z();
        sample3p_z();
    }

    for(float x = -500; x < 500; x+=fieldRes){
        for(float y = -500; y < 500; y+=fieldRes){
            for(float z = -500; z < 500; z+=fieldRes){
                // create a small red cube at each grid point
                cubes.emplace_back(vec3(x, y, z), vec4(1.0f, 0.0f, 0.0f, 0.05f), fieldRes);
                Cube &c = cubes[cubes.size() - 1];
                int count = 0;
                float half = c.size * 0.5f;
                float minX = c.position.x - half, maxX = c.position.x + half;
                float minY = c.position.y - half, maxY = c.position.y + half;
                float minZ = c.position.z - half, maxZ = c.position.z + half;

                std::unordered_map<uint32_t, int> colorCounts;
                for (const auto &p : particles) {
                    const vec3 &pos = p.position;
                    if (pos.x >= minX && pos.x < maxX &&
                        pos.y >= minY && pos.y < maxY &&
                        pos.z >= minZ && pos.z < maxZ) {
                        int r = int(p.color.r * 255.0f + 0.5f); r = std::max(0, std::min(255, r));
                        int g = int(p.color.g * 255.0f + 0.5f); g = std::max(0, std::min(255, g));
                        int b = int(p.color.b * 255.0f + 0.5f); b = std::max(0, std::min(255, b));
                        int a = int(p.color.a * 255.0f + 0.5f); a = std::max(0, std::min(255, a));
                        uint32_t key = (uint32_t(r) << 24) | (uint32_t(g) << 16) | (uint32_t(b) << 8) | uint32_t(a);
                        colorCounts[key] += 1;
                        ++count;
                    }
                }
                c.particleCount = count;

                // choose majority color (if any)
                if (count > 0) {
                    uint32_t bestKey = 0;
                    int bestCount = 0;
                    for (auto &kv : colorCounts) {
                        if (kv.second > bestCount) {
                            bestCount = kv.second;
                            bestKey = kv.first;
                        }
                    }
                    // unpack bytes back into normalized vec4
                    int br = (bestKey >> 24) & 0xFF;
                    int bg = (bestKey >> 16) & 0xFF;
                    int bb = (bestKey >> 8) & 0xFF;
                    int ba = (bestKey) & 0xFF;
                    c.color = vec4(br / 255.0f, bg / 255.0f, bb / 255.0f, ba / 255.0f);

                    // optionally scale alpha by density (keeps majority RGB)
                    // float alphaFromCount = std::min(1.0f, count / 500.0f);
                    // c.color.a = std::max(c.color.a, alphaFromCount);
                } else {
                    // empty cube: keep default low-alpha red
                    c.color = vec4(1.0f, 0.0f, 0.0f, 0.05f);
                }
                
                c.color.a = (count / 500.0f);
            }
        }
    }

    // -------- MAIN LOOP -------- //
    auto lastSampleTime = std::chrono::steady_clock::now();
    const std::chrono::milliseconds sampleInterval(100); // 0.1s

    while (!glfwWindowShouldClose(engine.window)) {
        engine.run();

        // ---- DRAW GRID ----
        grid.Draw(objectColorLoc);

        // ---- DRAW PARTICLES ----
        for (auto& p : particles) {
            p.Draw(objectColorLoc, modelLoc);
            glDrawArrays(GL_TRIANGLES, 0, p.vertices.size() / 3);
        }

        // ---- DRAW RED CUBES ----
        // draw transparent cubes back-to-front (simple approach: disable depth writes)
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glDepthMask(GL_FALSE); // don't write to depth buffer so blending works
        for (auto& c : cubes) {
            c.Draw(objectColorLoc, modelLoc);
        }
        glDepthMask(GL_TRUE);
        glDisable(GL_BLEND);

        glfwSwapBuffers(engine.window);
        glfwPollEvents();
    }

    // ---- CLEAN UP ----
    for (auto& p : particles) {
        glDeleteVertexArrays(1, &p.VAO);
        glDeleteBuffers(1, &p.VBO);
    }
    for (auto& c : cubes) {
        glDeleteVertexArrays(1, &c.VAO);
        glDeleteBuffers(1, &c.VBO);
    }
    glfwDestroyWindow(engine.window);
    glfwTerminate();
    return 0;
}