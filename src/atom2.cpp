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
const float a0 = 1;               //52.9f; // Bohr radius in pm
const float electron_r = 1.0f;
const double hbar  = 1;           //1.054571817e-34; // reduced Planck constant
const double m_e   = 1;           //9.10938356e-31;  // electron mass
const double zmSpeed   = 10.0;

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
    double zoomSpeed = zmSpeed;

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

        // Enable anti-aliased points for circular particles
        glEnable(GL_POINT_SMOOTH);

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
    float theta, phi;
    vec3 vel = vec3(0.0f);
    vec4 color = vec4(1.0f);
    Particle(vec3 p) : pos(p){}
    void drawParticle(GLint modelLoc, GLint objectColorLoc) {
        // Draw each particle
        glUniform4f(objectColorLoc, 0.5f, 1.0f, 1.0f, 20 / length(pos)); // Red color for particles
        glm::mat4 model = glm::mat4(1.0f);
        model = glm::translate(model, pos);
        model = glm::scale(model, glm::vec3(2.0f)); // Scale to make it visible
        glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model));

        // Make particle size larger when camera is closer for better visibility.
        float point_size = 100.0f / camera.radius;
        if (point_size < 1.0f) point_size = 1.0f;

        // Draw a simple point for the particle
        glPointSize(point_size);
        glBegin(GL_POINTS);
        glVertex3f(0.0f, 0.0f, 0.0f);
        glEnd();
    }
};

vector<Particle> particles;
random_device rd; mt19937 gen(rd()); uniform_real_distribution<float> dis(0.0f, 1.0f);

double sampleR(int n, int l, mt19937& gen) {
    const int N = 4096;
    //const double a0 = 1.0;
    const double rMax = 10.0 * n * n * a0;

    static vector<double> cdf;
    static bool built = false;

    if (!built) {
        cdf.resize(N);
        double dr = rMax / (N - 1);
        double sum = 0.0;

        for (int i = 0; i < N; ++i) {
            double r = i * dr;
            double rho = 2.0 * r / (n * a0);

            // Associated Laguerre L_{n-l-1}^{2l+1}(rho)
            int k = n - l - 1;
            int alpha = 2 * l + 1;

            double L = 1.0, Lm1 = 1.0 + alpha - rho;
            if (k == 1) L = Lm1;
            else if (k > 1) {
                double Lm2 = 1.0;
                for (int j = 2; j <= k; ++j) {
                    L = ((2*j - 1 + alpha - rho) * Lm1 -
                         (j - 1 + alpha) * Lm2) / j;
                    Lm2 = Lm1;
                    Lm1 = L;
                }
            }

            double norm = pow(2.0 / (n * a0), 3) * tgamma(n - l) / (2.0 * n * tgamma(n + l + 1));
            double R = sqrt(norm) * exp(-rho / 2.0) * pow(rho, l) * L;

            double pdf = r * r * R * R;
            sum += pdf;
            cdf[i] = sum;
        }

        for (double& v : cdf) v /= sum;
        built = true;
    }

    uniform_real_distribution<double> dis(0.0, 1.0);
    double u = dis(gen);

    int idx = lower_bound(cdf.begin(), cdf.end(), u) - cdf.begin();
    return idx * (rMax / (N - 1));
}
double sampleTheta(int l, int m, mt19937& gen) {
    const int N = 2048;
    static vector<double> cdf;
    static bool built = false;

    if (!built) {
        cdf.resize(N);
        double dtheta = M_PI / (N - 1);
        double sum = 0.0;

        for (int i = 0; i < N; ++i) {
            double theta = i * dtheta;
            double x = cos(theta);

            // Associated Legendre P_l^m(x)
            double Pmm = 1.0;
            if (m > 0) {
                double somx2 = sqrt((1.0 - x) * (1.0 + x));
                double fact = 1.0;
                for (int j = 1; j <= m; ++j) {
                    Pmm *= -fact * somx2;
                    fact += 2.0;
                }
            }

            double Plm;
            if (l == m) {
                Plm = Pmm;
            } else {
                double Pm1m = x * (2 * m + 1) * Pmm;
                if (l == m + 1) {
                    Plm = Pm1m;
                } else {
                    double Pll;
                    for (int ll = m + 2; ll <= l; ++ll) {
                        Pll = ((2 * ll - 1) * x * Pm1m -
                               (ll + m - 1) * Pmm) / (ll - m);
                        Pmm = Pm1m;
                        Pm1m = Pll;
                    }
                    Plm = Pm1m;
                }
            }

            double pdf = sin(theta) * Plm * Plm;
            sum += pdf;
            cdf[i] = sum;
        }

        for (double& v : cdf) v /= sum;
        built = true;
    }

    uniform_real_distribution<double> dis(0.0, 1.0);
    double u = dis(gen);

    int idx = lower_bound(cdf.begin(), cdf.end(), u) - cdf.begin();
    return idx * (M_PI / (N - 1));
}
float samplePhi(float n, float l, float m) {
    return 2.0f * M_PI * dis(gen);
}
vec3 calculateProbabilityFlow(Particle& p, int n, int l, int m) {
    double r = length(p.pos);   if (r < 1e-6) return vec3(0.0f);
    double theta = acos(p.pos.y / r); 
    double phi = atan2(p.pos.z, p.pos.x); 


    //Compute magnitude
    double sinTheta = sin(theta);  if (abs(sinTheta) < 1e-4) sinTheta = 1e-4;
    double v_mag = hbar * m / (m_e * r * sinTheta);

    //Convert to Cartesian
    double vx = -v_mag * sin(phi);
    double vy = 0.0; 
    double vz =  v_mag * cos(phi);

    return vec3((float)vx, (float)vy, (float)vz);
}

// ================= Main ================= //
int main () {
    setupCameraCallbacks(engine.window);
    GLint modelLoc = glGetUniformLocation(engine.shaderProgram, "model");
    GLint objectColorLoc = glGetUniformLocation(engine.shaderProgram, "objectColor");
    glUseProgram(engine.shaderProgram);

    float n = 6; float l = 4; float m = 1;

    // ------- CREATE PARTICLES -------
    for (int i = 0; i < 50000; ++i) {

        float r = sampleR(n, l, gen);
        float theta = sampleTheta(l, m, gen);
        float phi = samplePhi(n, l, m);

        vec3 pos = engine.sphericalToCartesian(r, theta, phi);
        particles.emplace_back(pos);
    }

    // ------------------ RENDERING LOOP ------------------
    float dt = 10.1f;
    while (!glfwWindowShouldClose(engine.window)) {
        engine.run();
        grid.Draw(objectColorLoc);

        for (Particle& p : particles) {

            // For a stationary state, Bohmian trajectories have constant r and theta.
            // We store the original r and theta before the update.
            double r = length(p.pos);
            if (r < 1e-6) { // Skip particles at the origin.
                p.drawParticle(modelLoc, objectColorLoc);
                continue;
            }
            double theta = acos(p.pos.y / r);
            p.vel = calculateProbabilityFlow(p, n, l, m);
            vec3 temp_pos = p.pos + p.vel * dt;
            double new_phi = atan2(temp_pos.z, temp_pos.x);
            p.pos = engine.sphericalToCartesian(r, theta, new_phi);

            if (p.pos.x > 0 || p.pos.y < 0)
            p.drawParticle(modelLoc, objectColorLoc);
        }

        glfwSwapBuffers(engine.window);
        glfwPollEvents();
    }
    
    glfwDestroyWindow(engine.window);
    glfwTerminate();
    return 0;
}




