#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <GL/glu.h>
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
const float a0 = 1;
float electron_r = 1.5f; // radius for spheres
const double hbar = 1;
const double m_e = 1;
const double zmSpeed = 10.0;

// --- Global quantum numbers ---
int n = 2, l = 1, m = 0, N = 100000;

// ================= Physics Sampling ================= //
struct Particle {
    vec3 pos;
    vec3 vel = vec3(0.0f);
    vec4 color;
    Particle(vec3 p, vec4 c = vec4(0.0f, 0.5f, 1.0f, 1.0f)) : pos(p), color(c){}
};
vector<Particle> particles;

// --- random devices ---
random_device rd; mt19937 gen(rd()); uniform_real_distribution<float> dis(0.0f, 1.0f);

// --- sample R ---                 <- uses CDF sampling
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
// --- sample Theta ---             <- uses CDF sampling
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
// --- sample Phi (uniform) ---     <- uses CDF sampling
float samplePhi(float n, float l, float m) {
    return 2.0f * M_PI * dis(gen);
}
// --- calculate prob current ---
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

// --- color map ---
vec4 heatmap_fire(float value) {
    // Ensure value is clamped between 0 and 1
    value = std::max(0.0f, std::min(1.0f, value));

    // Define color stops for the "Heat/Fire" pattern
    // Order: Black -> Dark Purple -> Red -> Orange -> Yellow -> White
    const int num_stops = 6;
    vec4 colors[num_stops] = {
        {0.0f, 0.0f, 0.0f, 1.0f}, // 0.0: Black
        {0.5f, 0.0f, 0.99f, 1.0f}, // 0.2: Dark Purple
        {0.8f, 0.0f, 0.0f, 1.0f}, // 0.4: Deep Red
        {1.0f, 0.5f, 0.0f, 1.0f}, // 0.6: Orange
        {1.0f, 1.0f, 0.0f, 1.0f}, // 0.8: Yellow
        {1.0f, 1.0f, 1.0f, 1.0f}  // 1.0: White
    };

    // Find which segment the value falls into
    float scaled_v = value * (num_stops - 1);
    int i = static_cast<int>(scaled_v);
    int next_i = std::min(i + 1, num_stops - 1);
    
    // Calculate how far we are between stop 'i' and 'next_i'
    float local_t = scaled_v - i;

    // Linearly interpolate between the two colors
    vec4 result;
    result.r = colors[i].r + local_t * (colors[next_i].r - colors[i].r);
    result.g = colors[i].g + local_t * (colors[next_i].g - colors[i].g);
    result.b = colors[i].b + local_t * (colors[next_i].b - colors[i].b);
    result.a = 1.0f; // Solid opacity

    return result;
}
vec4 inferno(double r, double theta, double phi, int n, int l, int m) {
    // --- radial part |R(r)|^2 ---
    double rho = 2.0 * r / (n * a0);

    int k = n - l - 1;
    int alpha = 2 * l + 1;

    double L = 1.0;
    if (k == 1) {
        L = 1.0 + alpha - rho;
    } else if (k > 1) {
        double Lm2 = 1.0;
        double Lm1 = 1.0 + alpha - rho;
        for (int j = 2; j <= k; ++j) {
            L = ((2*j - 1 + alpha - rho) * Lm1 -
                 (j - 1 + alpha) * Lm2) / j;
            Lm2 = Lm1;
            Lm1 = L;
        }
    }

    double norm = pow(2.0 / (n * a0), 3)
                * tgamma(n - l)
                / (2.0 * n * tgamma(n + l + 1));

    double R = sqrt(norm) * exp(-rho / 2.0) * pow(rho, l) * L;
    double radial = R * R;

    // --- angular part |P_l^m(cosθ)|^2 ---
    double x = cos(theta);

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
        double Pm1m = x * (2*m + 1) * Pmm;
        if (l == m + 1) {
            Plm = Pm1m;
        } else {
            for (int ll = m + 2; ll <= l; ++ll) {
                double Pll = ((2*ll - 1) * x * Pm1m -
                              (ll + m - 1) * Pmm) / (ll - m);
                Pmm = Pm1m;
                Pm1m = Pll;
            }
            Plm = Pm1m;
        }
    }

    double angular = Plm * Plm;

    double intensity = radial * angular;

    //cout << "intensity: " << intensity << endl;

    return heatmap_fire(intensity * 1.5 * pow(5, n)); // Scale for better color mapping
}


// ================= camera ================= //
struct Camera {
    vec3 target = vec3(0.0f, 0.0f, 0.0f);
    float radius = 50.0f;
    float azimuth = 0.0f;
    float elevation = M_PI / 2.0f;
    float orbitSpeed = 0.01f;
    float panSpeed = 0.01f;
    double zoomSpeed = zmSpeed;
    bool dragging = false;
    bool panning = false;
    double lastX = 0.0, lastY = 0.0;

    vec3 position() const {
        float clampedElevation = clamp(elevation, 0.01f, float(M_PI) - 0.01f);
        return vec3(
            radius * sin(clampedElevation) * cos(azimuth),
            radius * cos(clampedElevation),
            radius * sin(clampedElevation) * sin(azimuth)
        );
    }
    void update() {
        target = vec3(0.0f, 0.0f, 0.0f);
    }

    void processMouseMove(double x, double y) {
        float dx = float(x - lastX);
        float dy = float(y - lastY);
        if (dragging) {
            azimuth += dx * orbitSpeed;
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
                glfwGetCursorPos(win, &lastX, &lastY);
            } else if (action == GLFW_RELEASE) {
                dragging = false;
            }
        }
    }
    void processScroll(double xoffset, double yoffset) {
        radius -= yoffset * zoomSpeed;
        if (radius < 1.0f) radius = 1.0f;
        update();
    };
};
Camera camera;
vec3 sphericalToCartesian(float r, float theta, float phi){
        float x = r * sin(theta) * cos(phi);
        float y = r * cos(theta);
        float z = r * sin(theta) * sin(phi);
        return vec3(x, y, z);
    }
void generateParticles(int N) {
    particles.clear();
    for (int i = 0; i < N; ++i) {
        // --- get x, y, z, positions
        vec3 pos = sphericalToCartesian(
            sampleR(n, l, gen), 
            sampleTheta(l, m, gen), 
            samplePhi(n, l, m)
        );
        // --- color & add particle ---
        float r = length(pos);
        double theta = acos(pos.y / r);
        double phi = atan2(pos.z, pos.x);
        vec4 col = inferno(r, theta, phi, n, l, m) ;
        particles.emplace_back(pos, col);
    }
}

struct Engine {
    GLFWwindow* window;
    int WIDTH = 800;
    int HEIGHT = 600;

    // renders vars
    GLuint sphereVAO, sphereVBO;
    int sphereVertexCount;
    GLuint shaderProgram;
    GLint modelLoc, viewLoc, projLoc, colorLoc;

    // --- shaders ---
    const char* vertexShaderSource = R"glsl(
        #version 330 core
        layout(location=0) in vec3 aPos; uniform mat4 model; uniform mat4 view;
        uniform mat4 projection; out float lightIntensity;
        void main() { gl_Position = projection * view * model * vec4(aPos, 1.0);
            vec3 normal = normalize(aPos);
            vec3 lightDir = normalize(vec3(1.0, 1.0, 1.0));
            lightIntensity = max(dot(normal, lightDir), 0.5); // 0.2 is ambient light
        } )glsl";

    const char* fragmentShaderSource = R"glsl(
        #version 330 core
        in float lightIntensity; 
        out vec4 FragColor; 
        uniform vec4 objectColor;

        void main() { 
            // Increase the power to make the 'center-facing' spot tighter and brighter
            float glow = pow(lightIntensity, 2.0); 
            FragColor = vec4(objectColor.rgb , objectColor.a); 
        } )glsl";

    Engine() {
        if (!glfwInit()) exit(-1);
        window = glfwCreateWindow(800, 600, "Atom Prob-Flow", NULL, NULL);
        glfwMakeContextCurrent(window);
        glewInit();
        glEnable(GL_DEPTH_TEST);

        // Generate Sphere Vertices manually (like I did in the gravity sim)
        vector<float> vertices;
        float r = 0.05f; // Small sphere for particles
        int stacks = 10, sectors = 10;
        for(int i = 0; i <= stacks; ++i){
            float t1 = (float)i / stacks * M_PI;
            float t2 = (float)(i+1) / stacks * M_PI;
            for(int j = 0; j < sectors; ++j){
                float p1 = (float)j / sectors * 2 * M_PI;
                float p2 = (float)(j+1) / sectors * 2 * M_PI;
                auto getPos = [&](float t, float p) {
                    return vec3(r*sin(t)*cos(p), r*cos(t), r*sin(t)*sin(p));
                };
                vec3 v1 = getPos(t1, p1), v2 = getPos(t1, p2), v3 = getPos(t2, p1), v4 = getPos(t2, p2);
                vertices.insert(vertices.end(), {v1.x, v1.y, v1.z, v2.x, v2.y, v2.z, v3.x, v3.y, v3.z});
                vertices.insert(vertices.end(), {v2.x, v2.y, v2.z, v4.x, v4.y, v4.z, v3.x, v3.y, v3.z});
            }
        }
        sphereVertexCount = vertices.size() / 3;
        CreateVBOVAO(sphereVAO, sphereVBO, vertices);

        GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
        glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
        glCompileShader(vertexShader);

        GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
        glShaderSource(fragmentShader, 1, &fragmentShaderSource, NULL);
        glCompileShader(fragmentShader);

        shaderProgram = glCreateProgram();
        glAttachShader(shaderProgram, vertexShader);
        glAttachShader(shaderProgram, fragmentShader);
        glLinkProgram(shaderProgram);

        // Get uniform locations
        modelLoc = glGetUniformLocation(shaderProgram, "model");
        viewLoc  = glGetUniformLocation(shaderProgram, "view");
        projLoc  = glGetUniformLocation(shaderProgram, "projection");
        colorLoc = glGetUniformLocation(shaderProgram, "objectColor");
    }
    vec3 sphericalToCartesian(float r, float theta, float phi){
        float x = r * sin(theta) * cos(phi);
        float y = r * cos(theta);
        float z = r * sin(theta) * sin(phi);
        return vec3(x, y, z);
    }
    void CreateVBOVAO(GLuint& VAO, GLuint& VBO, const vector<float>& vertices) {
        glGenVertexArrays(1, &VAO);
        glGenBuffers(1, &VBO);
        glBindVertexArray(VAO);
        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_STATIC_DRAW);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);
    }
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
    void drawSpheres(vector<Particle>& particles) {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glUseProgram(shaderProgram); // Use our new shaded system

        mat4 projection = perspective(radians(45.0f), 800.0f/600.0f, 0.1f, 2000.0f);
        mat4 view = lookAt(camera.position(), camera.target, vec3(0, 1, 0)); 

        // Send view and projection to the shader
        glUniformMatrix4fv(viewLoc, 1, GL_FALSE, value_ptr(view));
        glUniformMatrix4fv(projLoc, 1, GL_FALSE, value_ptr(projection));

        glBindVertexArray(sphereVAO);

        for (auto& p : particles) {
            if (p.pos.x < 0 && p.pos.y > 0) continue;
            mat4 model = translate(mat4(1.0f), p.pos);
            model = scale(model, vec3(electron_r));
            glUniformMatrix4fv(modelLoc, 1, GL_FALSE, value_ptr(model));
            glUniform4f(colorLoc, p.color.r, p.color.g, p.color.b, p.color.a);
            
            glDrawArrays(GL_TRIANGLES, 0, sphereVertexCount);
        }
    }
    void setupCameraCallbacks() {
        glfwSetWindowUserPointer(window, &camera);
        glfwSetMouseButtonCallback(window, [](GLFWwindow* win, int button, int action, int mods) {
            ((Camera*)glfwGetWindowUserPointer(win))->processMouseButton(button, action, mods, win);
        });
        glfwSetCursorPosCallback(window, [](GLFWwindow* win, double x, double y) {
            ((Camera*)glfwGetWindowUserPointer(win))->processMouseMove(x, y);
        });
        glfwSetScrollCallback(window, [](GLFWwindow* win, double xoffset, double yoffset) {
            ((Camera*)glfwGetWindowUserPointer(win))->processScroll(xoffset, yoffset);
        });
        // Key callback: modify global quantum numbers
        glfwSetKeyCallback(window, [](GLFWwindow* win, int key, int scancode, int action, int mods) {
            if (!(action == GLFW_PRESS || action == GLFW_REPEAT)) return;

            if (key == GLFW_KEY_W) {
                n += 1;
                generateParticles(N);
            } else if (key == GLFW_KEY_S) {
                n -= 1;
                if (n < 1) n = 1;
                generateParticles(N);
            } else if (key == GLFW_KEY_E) {
                l += 1;
                generateParticles(N);
            } else if (key == GLFW_KEY_D) {
                l -= 1;
                if (l < 0) l = 0;
                generateParticles(N);
            } else if (key == GLFW_KEY_R) {
                m += 1;
                generateParticles(N);
            } else if (key == GLFW_KEY_F) {
                m -= 1;
                generateParticles(N);
            } else if (key == GLFW_KEY_T) {
                N +=100000;
                generateParticles(N);
            } else if (key == GLFW_KEY_G) {
                N -=100000;
                generateParticles(N);
            }

            // Clamp to valid ranges
            if (l > n - 1) l = n - 1;
            if (l < 0) l = 0;
            if (m > l) m = l;
            if (m < -l) m = -l;

            electron_r = float(n) / 3.0f;
            cout << "Quantum numbers updated: n=" << n << " l=" << l << " m=" << m << " N=" << N << "\n";
        });
    }
};
Engine engine;

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

// ================= Main Loop ================= //
int main () {
    GLint modelLoc = glGetUniformLocation(engine.shaderProgram, "model");
    GLint objectColorLoc = glGetUniformLocation(engine.shaderProgram, "objectColor");
    glUseProgram(engine.shaderProgram);
    engine.setupCameraCallbacks();

    // --- scale r for bigger orbitals ---
    electron_r = float(n) / 3.0f;

    // --- Sample particles ---
    generateParticles(250000);

    float dt = 0.5f;
    cout << "Starting simulation..." << endl;
    while (!glfwWindowShouldClose(engine.window)) {
        grid.Draw(objectColorLoc);

        // ------ Update Probability current ------
        for (Particle& p : particles) {
            double r = length(p.pos);
            if (r > 1e-6) {
                double theta = acos(p.pos.y / r);
                p.vel = calculateProbabilityFlow(p, n, l, m);
                vec3 temp_pos = p.pos + p.vel * dt;
                double new_phi = atan2(temp_pos.z, temp_pos.x);
                p.pos = engine.sphericalToCartesian(r, theta, new_phi);
            }
        }
        // ------ Draw Particles ------
        engine.drawSpheres(particles);

        glfwSwapBuffers(engine.window);
        glfwPollEvents();
    }

    // --- close ---
    glfwDestroyWindow(engine.window);
    glfwTerminate();
    return 0;
}