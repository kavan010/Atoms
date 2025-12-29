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
#include <fstream>
#include <complex>
#include <random>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
using namespace glm;
using namespace std;

// ================= Constants ================= //
const float c = 299792458.0f / 100000000.0f;    // speed of light in m/s
const float a0 = 52.9f; // Bohr radius in pm
const float electron_r = 1.0f;
const double hbar = 1.054571817e-34; // reduced Planck constant
const double m_e   = 9.10938356e-31;  // electron mass

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
    double zoomSpeed = 20.0f;

    bool dragging = false;
    bool panning = false;
    bool moving = false; // For compute shader optimization
    double lastX = 0.0, lastY = 0.0;

    // Calculate camera position in world space
    vec3 position() const {
        float clampedElevation = glm::clamp(elevation, 0.01f, float(M_PI) - 0.01f);
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
            elevation = glm::clamp(elevation, 0.01f, float(M_PI) - 0.01f);
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

// ================= Physics ================= //
struct Physics {
    // Factorial helper
    double factorial(int n) {
        if (n <= 1) return 1.0;
        double res = 1.0;
        for (int i = 2; i <= n; ++i) res *= i;
        return res;
    }
    // Associated Laguerre Polynomial L_n^alpha(x)
    double laguerre(int n, int alpha, double x) {
        if (n == 0) return 1.0;
        if (n == 1) return 1.0 + alpha - x;
        
        double L_k_minus_1 = 1.0 + alpha - x;
        double L_k_minus_2 = 1.0;
        double L_k = 0.0;

        for (int k = 1; k < n; ++k) {
            L_k = ((2 * k + 1 + alpha - x) * L_k_minus_1 - (k + alpha) * L_k_minus_2) / (k + 1);
            L_k_minus_2 = L_k_minus_1;
            L_k_minus_1 = L_k;
        }
        return L_k;
    }
    // Associated Legendre Polynomial P_l^m(x)
    double legendre(int l, int m, double x) {
        // Condon-Shortley phase is often included in Spherical Harmonics, 
        // but for probability density |Y|^2 it doesn't matter.
        // Here we use a standard recurrence.
        
        double pmm = 1.0;
        if (m > 0) {
            double somx2 = sqrt((1.0 - x) * (1.0 + x));
            double fact = 1.0;
            for (int i = 1; i <= m; ++i) {
                pmm *= -fact * somx2;
                fact += 2.0;
            }
        }
        if (l == m) return pmm;

        double pmmp1 = x * (2 * m + 1) * pmm;
        if (l == m + 1) return pmmp1;

        double pll = 0.0;
        for (int ll = m + 2; ll <= l; ++ll) {
            pll = (x * (2 * ll - 1) * pmmp1 - (ll + m - 1) * pmm) / (ll - m);
            pmm = pmmp1;
            pmmp1 = pll;
        }
        return pll;
    }
    // Hydrogen Radial Function R_nl(r)
    double radial_R(int n, int l, double r) {
        double rho = 2.0 * r / (n * a0);
        double prefactor = sqrt(pow(2.0 / (n * a0), 3) * factorial(n - l - 1) / (2.0 * n * factorial(n + l)));
        return prefactor * exp(-rho / 2.0) * pow(rho, l) * laguerre(n - l - 1, 2 * l + 1, rho);
    }
    // Probability Density Function |Psi|^2
    double probability_density(int n, int l, int m, const vec3& pos) {
        double r = length(pos);
        if (r < 1e-6) return 0.0;

        // Convert Cartesian to Spherical
        double theta = acos(pos.z / r); // Elevation [0, PI]
        // phi is not needed for magnitude of stationary states as |e^im*phi| = 1
        
        double R = radial_R(n, l, r);
        
        // Normalization for Spherical Harmonics (part of it)
        double Y_norm = sqrt(((2 * l + 1) * factorial(l - abs(m))) / (4 * M_PI * factorial(l + abs(m))));
        double P = legendre(l, abs(m), cos(theta));
        
        double psi_mag = R * Y_norm * P;
        return psi_mag * psi_mag;
    }
    // Bohmian Velocity Calculation
    vec3 getBohmianVelocity(int n, int l, int m, const vec3& pos) {
        // For a single eigenstate psi_{nlm}, the flow is purely azimuthal.
        // v = (hbar / m_e) * (m / (r * sin(theta))) in phi_hat direction
        
        double r = length(pos);
        if (r < 1e-4) return vec3(0.0);

        // r_xy is r * sin(theta) (distance from Z axis)
        double r_xy = sqrt(pos.x * pos.x + pos.y * pos.y); 
        if (r_xy < 1e-4) return vec3(0.0); // Singularity at poles

        // Velocity magnitude
        // We scale this down significantly for visual stability
        double speed = (m) / r_xy; // Simplified unit-less speed relative to shape
        
        // Phi unit vector: (-y, x, 0) / r_xy
        vec3 phi_hat = vec3(-pos.y / r_xy, pos.x / r_xy, 0.0f);

        return phi_hat * (float)speed;
    }
};
Physics phy;

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
        glUniform4f(objectColorLoc, 1.0f, 1.0f, 1.0f, 0.5f);
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

struct Particle {
    vec3 pos;
    vec3 vel = vec3(0.0f);
    vec4 color = vec4(0.5f, 1.0f, 1.0f, 0.9f);
    Particle(vec3 p) : pos(p){}
    void drawParticle(GLint modelLoc, GLint objectColorLoc) {
        // Draw each particle
        glUniform4f(objectColorLoc, 0.5f, 1.0f, 1.0f, 0.9f); // Red color for particles
        glm::mat4 model = glm::mat4(1.0f);
        model = glm::translate(model, pos);
        model = glm::scale(model, glm::vec3(2.0f)); // Scale to make it visible
        glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model));

        // Draw a simple point for the particle
        glPointSize(electron_r);
        glBegin(GL_POINTS);
        glVertex3f(0.0f, 0.0f, 0.0f);
        glEnd();
    }
};

vector<Particle> particles;
vector<Particle> GenerateParticles(int n, int l, int m, int count) {
    vector<Particle> particles;
    particles.reserve(count);

    // Initial random walker position
    vec3 walker = vec3(n*a0 + 1.0, 0, 0); 
    double currentProb = phy.probability_density(n, l, m, walker);

    default_random_engine generator(time(0));
    normal_distribution<double> stepDist(0.0, 1.5); // Step size for walker
    uniform_real_distribution<double> acceptDist(0.0, 1.0);

    // Warmup (burn-in) to find the cloud
    for(int i=0; i<5000; i++) {
        vec3 proposal = walker + vec3(stepDist(generator), stepDist(generator), stepDist(generator));
        double nextProb = phy.probability_density(n, l, m, proposal);
        if (nextProb > 0 && acceptDist(generator) < (nextProb / currentProb)) {
            walker = proposal;
            currentProb = nextProb;
        }
    }

    // Actual Sampling
    int accepted = 0;
    while(accepted < count) {
        // Decorrelation steps (move walker a few times between samples to avoid clumping)
        for(int k=0; k<10; k++) {
            vec3 proposal = walker + vec3(stepDist(generator), stepDist(generator), stepDist(generator));
            double nextProb = phy.probability_density(n, l, m, proposal);
            // Metropolis acceptance criterion
            if (nextProb >= currentProb || acceptDist(generator) < (nextProb / currentProb)) {
                walker = proposal;
                currentProb = nextProb;
            }
        }
        
        // Color based on phase (simplified for real/im split or just depth)
        float r_dist = length(walker);
        float brightness = 0.2f + 0.8f * (1.0f - exp(-r_dist/10.0f));
        vec4 col = vec4(0.2f, 0.6f, 1.0f, 0.6f); 

        particles.push_back(Particle(walker));
        accepted++;
        
        if (accepted % 1000 == 0) cout << "Sampled " << accepted << " particles..." << endl;
    }
    return particles;
}

// ================= Main ================= //
int main () {
    setupCameraCallbacks(engine.window);
    GLint modelLoc = glGetUniformLocation(engine.shaderProgram, "model");
    GLint objectColorLoc = glGetUniformLocation(engine.shaderProgram, "objectColor");
    glUseProgram(engine.shaderProgram);

    float n = 6; float l = 4; float m = 1;

    // ------- CREATE PARTICLES -------
    particles = GenerateParticles(n, l, m, 10000);

    // ------------------ RENDERING LOOP ------------------
    float dt = 0.1f;
    while (!glfwWindowShouldClose(engine.window)) {
        engine.run();
        grid.Draw(objectColorLoc);

        for (Particle& p : particles) {

            // Standard RK4
            vec3 k1 = phy.getBohmianVelocity(n, l, m, p.pos);
            vec3 k2 = phy.getBohmianVelocity(n, l, m, p.pos + k1 * (dt * 0.5f));
            vec3 k3 = phy.getBohmianVelocity(n, l, m, p.pos + k2 * (dt * 0.5f));
            vec3 k4 = phy.getBohmianVelocity(n, l, m, p.pos + k3 * dt);

            p.pos += (k1 + 2.0f*k2 + 2.0f*k3 + k4) * (dt / 6.0f) ;

            p.drawParticle(modelLoc, objectColorLoc);
        }

        glfwSwapBuffers(engine.window);
        glfwPollEvents();
    }
    
    glfwDestroyWindow(engine.window);
    glfwTerminate();
    return 0;
}