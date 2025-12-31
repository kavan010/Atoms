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
float electron_r = 1.5f; // Will be used as radius for ray-traced spheres
const double hbar = 1;
const double m_e = 1;
const double zmSpeed = 10.0;

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

// --- sample R ---
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
// --- sample Theta ---
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
// --- sample Phi (uniform) ---
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

    return heatmap_fire(intensity * 1500); // Scale for better color mapping
}


// ================= Raytracer ================= //
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

struct Engine {
    GLFWwindow* window;
    int WIDTH = 800;
    int HEIGHT = 600;

    // Raytracing vals
    GLuint sphereVAO, sphereVBO;
    int sphereVertexCount;
    GLuint shaderProgram;
    GLint modelLoc, viewLoc, projLoc, colorLoc;

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

        // Generate Sphere Vertices manually (like you did in the gravity sim)
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
    }
};
Engine engine;



// ================= Main Loop ================= //
int main () {
    engine.setupCameraCallbacks();

    // --- Quantum numbers setup ---
    float n = 5; float l = 4; float m = 1;
    electron_r = n / 3 ;

    // --- Sample particles ---
    for (int i = 0; i < 500000; ++i) {
        // --- get x, y, z, positions
        vec3 pos = engine.sphericalToCartesian(
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

    // Inside main(), before the while loop:
    GLfloat spec[] = { 1.0f, 1.0f, 1.0f, 1.0f };
    GLfloat shininess[] = { 50.0f };
    glMaterialfv(GL_FRONT, GL_SPECULAR, spec);
    glMaterialfv(GL_FRONT, GL_SHININESS, shininess);

    float dt = 0.5f;
    cout << "Starting simulation..." << endl;
    while (!glfwWindowShouldClose(engine.window)) {

        // ------ Draw Particles ------
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
        engine.drawSpheres(particles);

        glfwSwapBuffers(engine.window);
        glfwPollEvents();
    }

    // --- Cleanup ---
    
    glfwDestroyWindow(engine.window);
    glfwTerminate();
    return 0;
}









// #include <GL/glew.h>
// #include <GLFW/glfw3.h>
// #include <glm/glm.hpp>
// #include <glm/gtc/matrix_transform.hpp>
// #include <glm/gtc/type_ptr.hpp>
// #include <vector>
// #include <iostream>

// const char* vertexShaderSource = R"glsl(
// #version 330 core
// layout(location=0) in vec3 aPos;
// uniform mat4 model;
// uniform mat4 view;
// uniform mat4 projection;
// out float lightIntensity;
// void main() {
//     gl_Position = projection * view * model * vec4(aPos, 1.0);
//     vec3 worldPos = (model * vec4(aPos, 1.0)).xyz;
//     vec3 normal = normalize(aPos);
//     vec3 dirToCenter = normalize(-worldPos);
//     lightIntensity = max(dot(normal, dirToCenter), 0.3);})glsl";

// const char* fragmentShaderSource = R"glsl(
// #version 330 core
// in float lightIntensity;
// out vec4 FragColor;
// uniform vec4 objectColor;
// uniform bool isGrid; // Add this uniform
// uniform bool GLOW;
// void main() {
//     if (isGrid) {
//         // If it's the grid, use the original color without lighting
//         FragColor = objectColor;
//     } else if(GLOW){
//         FragColor = vec4(objectColor.rgb * 10000000, objectColor.a);
//     }else {
//         // If it's an object, apply the lighting effect
//         float fade = smoothstep(0.0, 10.0, lightIntensity*10);
//         FragColor = vec4(objectColor.rgb * fade, objectColor.a);
//     }})glsl";

// bool running = true;
// bool pause = true;
// float deltaTime = 0.0;
// float lastFrame = 0.0;

// const double G = 6.6743e-11; // m^3 kg^-1 s^-2
// const float c = 299792458.0;
// float initMass = float(pow(10, 23));
// float sizeRatio = 30000.0f;
// GLFWwindow* StartGLU();
// GLuint CreateShaderProgram(const char* vertexSource, const char* fragmentSource);
// void CreateVBOVAO(GLuint& VAO, GLuint& VBO, const float* vertices, size_t vertexCount);
// glm::vec3 sphericalToCartesian(float r, float theta, float phi);

// struct Camera{
//     // Camera vals
//     float distance;
//     float yaw;
//     float pitch;
//     glm::vec3 position;
//     glm::vec3 target;
//     glm::vec3 up;

//     // mouse vals
//     bool middleMousePressed = false;
//     double lastX = 0.0, lastY = 0.0;
//     float orbitSpeed = 0.4f;
//     float zoomSpeed = 1500.0f;

//     float fov = 60.0f; 
//     double lastMovementTime = 0.0;

//     Camera(glm::vec3 t = glm::vec3(0.0f, 0.0f, -19.0f), float dist = 5.0f, float yawVal = -90.0f, float pitchVal = 0.0f, float fovVal = 90.0f)
//         : target(t), distance(dist), yaw(yawVal), pitch(pitchVal), fov(fovVal) {
//         up = glm::vec3(0, 1, 0);
//     }

//     void updatePosition() {
//         float radYaw = glm::radians(yaw);
//         float radPitch = glm::radians(pitch);
//         position.x = target.x + distance * cos(radPitch) * cos(radYaw);
//         position.y = target.y + distance * sin(radPitch);
//         position.z = target.z + distance * cos(radPitch) * sin(radYaw);
//     }

//     // member funcs
//     void handleMidMouse(int button, int action, int mods, GLFWwindow* window) {
//         if (button == GLFW_MOUSE_BUTTON_MIDDLE) {
//             if (action == GLFW_PRESS) {
//                 middleMousePressed = true;
//                 glfwGetCursorPos(window, &lastX, &lastY);
//             } else if (action == GLFW_RELEASE) {
//                 middleMousePressed = false;
//             }
//         }
//     }
//     void handleCursorPos(double xpos, double ypos, GLFWwindow* window) {
//         if (!middleMousePressed)
//             return;

//         double deltaX = xpos - lastX;
//         double deltaY = -(ypos - lastY);

//         // If shift, pan camera
//         if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS ||
//             glfwGetKey(window, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS) {
//             glm::vec3 forward = normalize(target - position);
//             glm::vec3 right = normalize(cross(forward, up));
//             glm::vec3 camUp = cross(right, forward);
//             float panSpeed = 0.005f * distance;
//             target += -right * (float)deltaX * panSpeed + camUp * (float)deltaY * panSpeed;
//         }
//         // orbit camera.
//         else {
//             yaw   += (float)deltaX * orbitSpeed;
//             pitch += (float)deltaY * orbitSpeed;
//             if (pitch > 89.0f)  pitch = 89.0f;
//             if (pitch < -89.0f) pitch = -89.0f;
//         }
//         updatePosition();
//         lastX = xpos;
//         lastY = ypos;
//         lastMovementTime = glfwGetTime();
//     }
//     void handleScroll(double xoffset, double yoffset, GLFWwindow* window) {

//         distance -= (float)yoffset * zoomSpeed;
        
//         updatePosition();
//         lastMovementTime = glfwGetTime();
//     }
//     // keys
//     void handleKeys(GLFWwindow* window, int key, int scancode, int action, int mods) {

//     }

//     // boring camera functions D:
//     static void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods) {
//         Camera* cam = static_cast<Camera*>(glfwGetWindowUserPointer(window));
//         cam->handleMidMouse(button, action, mods, window);
//     }
//     static void cursorPositionCallback(GLFWwindow* window, double xpos, double ypos) {
//         Camera* cam = static_cast<Camera*>(glfwGetWindowUserPointer(window));
//         cam->handleCursorPos(xpos, ypos, window);
//     }
//     static void scrollCallback(GLFWwindow* window, double xoffset, double yoffset) {
//         Camera* cam = static_cast<Camera*>(glfwGetWindowUserPointer(window));
//         cam->handleScroll(xoffset, yoffset, window);
//     }
//     static void KeyCallback(GLFWwindow* window, int key, int scancode, int action, int mods) {
//         Camera* cam = static_cast<Camera*>(glfwGetWindowUserPointer(window));
//         cam->handleKeys(window, key, scancode, action, mods);
//     }

//     void registerCallbacks(GLFWwindow* window) {
//         glfwSetWindowUserPointer(window, this);
//         glfwSetMouseButtonCallback(window, Camera::mouseButtonCallback);
//         glfwSetKeyCallback(window, Camera::KeyCallback);
//         glfwSetCursorPosCallback(window, Camera::cursorPositionCallback);
//         glfwSetScrollCallback(window, Camera::scrollCallback);
//     }
// };
// class Object {
//     public:
//         GLuint VAO, VBO;
//         glm::vec3 position = glm::vec3(400, 300, 0);
//         glm::vec3 velocity = glm::vec3(0, 0, 0);
//         size_t vertexCount;
//         glm::vec4 color = glm::vec4(1.0f, 0.0f, 0.0f, 1.0f);

//         bool Initalizing = false;
//         bool Launched = false;
//         bool target = false;

//         float mass;
//         float density;  // kg / m^3  HYDROGEN
//         float radius;

//         glm::vec3 LastPos = position;
//         bool glow;

//         Object(glm::vec3 initPosition, glm::vec3 initVelocity, float mass, float density = 3344, glm::vec4 color = glm::vec4(1.0f, 0.0f, 0.0f, 1.0f), bool Glow = false) {   
//             this->position = initPosition;
//             this->velocity = initVelocity;
//             this->mass = mass;
//             this->density = density;
//             this->radius = pow(((3 * this->mass/this->density)/(4 * 3.14159265359)), (1.0f/3.0f)) / sizeRatio;
//             this->color = color;
//             this->glow = Glow;
            

//             // Generate vertices (centered at origin)
//             std::vector<float> vertices = Draw();
//             vertexCount = vertices.size();

//             CreateVBOVAO(VAO, VBO, vertices.data(), vertexCount);
//         }

//         std::vector<float> Draw() {
//             std::vector<float> vertices;
//             int stacks = 25;
//             int sectors = 25;
            
//             for(float i = 0.0f; i <= stacks; ++i){
//                 float theta1 = (i / stacks) * glm::pi<float>();
//                 float theta2 = (i+1) / stacks * glm::pi<float>();
//                 for (float j = 0.0f; j < sectors; ++j){
//                     float phi1 = j / sectors * 2 * glm::pi<float>();
//                     float phi2 = (j+1) / sectors * 2 * glm::pi<float>();
//                     glm::vec3 v1 = sphericalToCartesian(this->radius, theta1, phi1);
//                     glm::vec3 v2 = sphericalToCartesian(this->radius, theta1, phi2);
//                     glm::vec3 v3 = sphericalToCartesian(this->radius, theta2, phi1);
//                     glm::vec3 v4 = sphericalToCartesian(this->radius, theta2, phi2);

//                     // Triangle 1: v1-v2-v3
//                     vertices.insert(vertices.end(), {v1.x, v1.y, v1.z}); //      /|
//                     vertices.insert(vertices.end(), {v2.x, v2.y, v2.z}); //     / |
//                     vertices.insert(vertices.end(), {v3.x, v3.y, v3.z}); //    /__|
                    
//                     // Triangle 2: v2-v4-v3
//                     vertices.insert(vertices.end(), {v2.x, v2.y, v2.z});
//                     vertices.insert(vertices.end(), {v4.x, v4.y, v4.z});
//                     vertices.insert(vertices.end(), {v3.x, v3.y, v3.z});
//                 }   
//             }
//             return vertices;
//         }
        
//         void UpdatePos(){
//             this->position[0] += this->velocity[0] / 94;
//             this->position[1] += this->velocity[1] / 94;
//             this->position[2] += this->velocity[2] / 94;
//             this->radius = pow(((3 * this->mass/this->density)/(4 * 3.14159265359)), (1.0f/3.0f)) / sizeRatio;
//         }
//         void UpdateVertices() {
//             // Generate new vertices with current radius
//             std::vector<float> vertices = Draw();
            
//             // Update the VBO with new vertex data
//             glBindBuffer(GL_ARRAY_BUFFER, VBO);
//             glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_STATIC_DRAW);
//         }
//         glm::vec3 GetPos() const {
//             return this->position;
//         }
//         void accelerate(float x, float y, float z){
//             this->velocity[0] += x / 96;
//             this->velocity[1] += y / 96;
//             this->velocity[2] += z / 96;
//         }
//         float CheckCollision(const Object& other) {
//             float dx = other.position[0] - this->position[0];
//             float dy = other.position[1] - this->position[1];
//             float dz = other.position[2] - this->position[2];
//             float distance = std::pow(dx*dx + dy*dy + dz*dz, (1.0f/2.0f));
//             // if (other.radius + this->radius > distance){
//             //     return -0.2f;
//             // }
//             return 1.0f;
//         }
// };

// std::vector<float> CreateGridVertices(float size, int divisions, const std::vector<Object>& objs);
// std::vector<float> UpdateGridVertices(std::vector<float> vertices, const std::vector<Object>& objs, float halfSize, float originalY);
// void DrawGrid(GLuint shaderProgram, GLuint gridVAO, size_t vertexCount);
// void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods);
// void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods);
// void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);

// std::vector<Object> objs = {};

// GLuint gridVAO, gridVBO;

// // Instantiate the Camera class
// Camera camera(glm::vec3(0.0f, 5000.0f, 5000.0f), 10000.0f, -90.0f, 0.0f, 45.0f);

// int main() {
//     GLFWwindow* window = StartGLU();
//     GLuint shaderProgram = CreateShaderProgram(vertexShaderSource, fragmentShaderSource);

//     GLint modelLoc = glGetUniformLocation(shaderProgram, "model");
//     GLint objectColorLoc = glGetUniformLocation(shaderProgram, "objectColor");
//     glUseProgram(shaderProgram);

//     // Register Camera Callbacks
//     camera.registerCallbacks(window);

//     // Projection matrix
//     glm::mat4 projection = glm::perspective(glm::radians(camera.fov), 800.0f / 600.0f, 0.1f, 750000.0f);
//     GLint projectionLoc = glGetUniformLocation(shaderProgram, "projection");
//     glUniformMatrix4fv(projectionLoc, 1, GL_FALSE, glm::value_ptr(projection));

//     objs = {
//         Object(glm::vec3(0, 0, 0), glm::vec3(0, 0, 0), 191000000000000000000000000000.0f, 208000000000.0f, glm::vec4(1.0f, 0.929f, 0.176f, 1.0f), true),
//         Object(glm::vec3(0, 0, 10000), glm::vec3(0, 0, 0), 191000000000000000000000000000.0f, 208000000.0f, glm::vec4(1.0f, 0.929f, 0.176f, 1.0f), true),
//     };

//     float size = 40000.0f;
//     int divisions = 50;
//     float step = size / divisions;
//     float halfSize = size / 2.0f;
//     float originalY = -halfSize * 0.3f + 3 * step;

//     std::vector<float> gridVertices = CreateGridVertices(size, divisions, objs);
//     CreateVBOVAO(gridVAO, gridVBO, gridVertices.data(), gridVertices.size());

//     while (!glfwWindowShouldClose(window) && running == true) {
//         float currentFrame = glfwGetTime();
//         deltaTime = currentFrame - lastFrame;
//         lastFrame = currentFrame;

//         glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

//         glfwSetKeyCallback(window, camera.KeyCallback);
//         glfwSetMouseButtonCallback(window, camera.mouseButtonCallback);

//         // Update Camera View Matrix
//         camera.updatePosition();
//         glm::mat4 view = glm::lookAt(camera.position, camera.target, camera.up);
//         GLint viewLoc = glGetUniformLocation(shaderProgram, "view");
//         glUniformMatrix4fv(viewLoc, 1, GL_FALSE, glm::value_ptr(view));

//         // Update objects initializing
//         if (!objs.empty() && objs.back().Initalizing) {
//             if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS) {
//                 objs.back().mass *= 1.0 + 5.0 * deltaTime;
//                 std::cout << "radius: " << objs.back().radius << std::endl;

//                 // Update vertex data
//                 objs.back().UpdateVertices();
//             }
//         }

//         // Draw the grid
//         glUseProgram(shaderProgram);
//         glUniform4f(objectColorLoc, 1.0f, 1.0f, 1.0f, 0.25f);
//         glUniform1i(glGetUniformLocation(shaderProgram, "isGrid"), 1);
//         glUniform1i(glGetUniformLocation(shaderProgram, "GLOW"), 0);
//         gridVertices = UpdateGridVertices(gridVertices, objs, halfSize, originalY);
//         glBindBuffer(GL_ARRAY_BUFFER, gridVBO);
//         glBufferData(GL_ARRAY_BUFFER, gridVertices.size() * sizeof(float), gridVertices.data(), GL_DYNAMIC_DRAW);
//         DrawGrid(shaderProgram, gridVAO, gridVertices.size());

//         // Draw the objects
//         for (auto& obj : objs) {
//             glUniform4f(objectColorLoc, obj.color.r, obj.color.g, obj.color.b, obj.color.a);

//             for (auto& obj2 : objs) {
//                 if (&obj2 != &obj && !obj.Initalizing && !obj2.Initalizing && !obj.Launched && !obj2.Launched) {
//                     float dx = obj2.GetPos()[0] - obj.GetPos()[0];
//                     float dy = obj2.GetPos()[1] - obj.GetPos()[1];
//                     float dz = obj2.GetPos()[2] - obj.GetPos()[2];
//                     float distance = sqrt(dx * dx + dy * dy + dz * dz);

//                     if (distance > 0) {
//                         std::vector<float> direction = {dx / distance, dy / distance, dz / distance};
//                         distance *= 1000;
//                         double Gforce = (G * obj.mass * obj2.mass) / (distance * distance);

//                         float acc1 = Gforce / obj.mass;
//                         std::vector<float> acc = {direction[0] * acc1, direction[1] * acc1, direction[2] * acc1};
//                         if (!pause) {
//                             obj.accelerate(acc[0], acc[1], acc[2]);
//                         }

//                         // Collision
//                         obj.velocity *= obj.CheckCollision(obj2);
//                     }
//                 }
//             }
//             if (obj.Initalizing) {
//                 obj.radius = pow(((3 * obj.mass / obj.density) / (4 * 3.14159265359)), (1.0f / 3.0f)) / sizeRatio;
//                 obj.UpdateVertices();
//                 obj.glow = true;
//             }

//             if (obj.Launched) {
//                 obj.Launched = false;
//             }

//             // Update positions
//             if (!pause) {
//                 obj.UpdatePos();
//             }

//             glm::mat4 model = glm::mat4(1.0f);
//             model = glm::translate(model, obj.position); // Apply position here
//             glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model));
//             glUniform1i(glGetUniformLocation(shaderProgram, "isGrid"), 0);
//             if (obj.glow) {
//                 glUniform1i(glGetUniformLocation(shaderProgram, "GLOW"), 1);
//             } else {
//                 glUniform1i(glGetUniformLocation(shaderProgram, "GLOW"), 0);
//             }

//             glBindVertexArray(obj.VAO);
//             glDrawArrays(GL_TRIANGLES, 0, obj.vertexCount / 3);
//         }

//         glfwSwapBuffers(window);
//         glfwPollEvents();
//     }

//     for (auto& obj : objs) {
//         glDeleteVertexArrays(1, &obj.VAO);
//         glDeleteBuffers(1, &obj.VBO);
//     }

//     glDeleteVertexArrays(1, &gridVAO);
//     glDeleteBuffers(1, &gridVBO);

//     glDeleteProgram(shaderProgram);
//     glfwTerminate();

//     return 0;
// }

// GLFWwindow* StartGLU() {
//     if (!glfwInit()) {
//         std::cout << "Failed to initialize GLFW, panic" << std::endl;
//         return nullptr;
//     }
//     GLFWwindow* window = glfwCreateWindow(800, 600, "3D_TEST", NULL, NULL);
//     if (!window) {
//         std::cerr << "Failed to create GLFW window." << std::endl;
//         glfwTerminate();
//         return nullptr;
//     }
//     glfwMakeContextCurrent(window);

//     glewExperimental = GL_TRUE;
//     if (glewInit() != GLEW_OK) {
//         std::cerr << "Failed to initialize GLEW." << std::endl;
//         glfwTerminate();
//         return nullptr;
//     }

//     glEnable(GL_DEPTH_TEST);
//     glViewport(0, 0, 800, 600);
//     glEnable(GL_BLEND);
//     glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); // Standard blending for transparency

//     return window;
// }

// GLuint CreateShaderProgram(const char* vertexSource, const char* fragmentSource) {
//     // Vertex shader
//     GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
//     glShaderSource(vertexShader, 1, &vertexSource, nullptr);
//     glCompileShader(vertexShader);

//     GLint success;
//     glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success);
//     if (!success) {
//         char infoLog[512];
//         glGetShaderInfoLog(vertexShader, 512, nullptr, infoLog);
//         std::cerr << "Vertex shader compilation failed: " << infoLog << std::endl;
//     }

//     // Fragment shader
//     GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
//     glShaderSource(fragmentShader, 1, &fragmentSource, nullptr);
//     glCompileShader(fragmentShader);

//     glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &success);
//     if (!success) {
//         char infoLog[512];
//         glGetShaderInfoLog(fragmentShader, 512, nullptr, infoLog);
//         std::cerr << "Fragment shader compilation failed: " << infoLog << std::endl;
//     }

//     // Shader program
//     GLuint shaderProgram = glCreateProgram();
//     glAttachShader(shaderProgram, vertexShader);
//     glAttachShader(shaderProgram, fragmentShader);
//     glLinkProgram(shaderProgram);

//     glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);
//     if (!success) {
//         char infoLog[512];
//         glGetProgramInfoLog(shaderProgram, 512, nullptr, infoLog);
//         std::cerr << "Shader program linking failed: " << infoLog << std::endl;
//     }

//     glDeleteShader(vertexShader);
//     glDeleteShader(fragmentShader);

//     return shaderProgram;
// }
// void CreateVBOVAO(GLuint& VAO, GLuint& VBO, const float* vertices, size_t vertexCount) {
//     glGenVertexArrays(1, &VAO);
//     glGenBuffers(1, &VBO);

//     glBindVertexArray(VAO);
//     glBindBuffer(GL_ARRAY_BUFFER, VBO);
//     glBufferData(GL_ARRAY_BUFFER, vertexCount * sizeof(float), vertices, GL_STATIC_DRAW);

//     glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
//     glEnableVertexAttribArray(0);
//     glBindVertexArray(0);
// }

// glm::vec3 sphericalToCartesian(float r, float theta, float phi){
//     float x = r * sin(theta) * cos(phi);
//     float y = r * cos(theta);
//     float z = r * sin(theta) * sin(phi);
//     return glm::vec3(x, y, z);
// };
// void DrawGrid(GLuint shaderProgram, GLuint gridVAO, size_t vertexCount) {
//     glUseProgram(shaderProgram);
//     glm::mat4 model = glm::mat4(1.0f); // Identity matrix for the grid
//     GLint modelLoc = glGetUniformLocation(shaderProgram, "model");
//     glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model));

//     glBindVertexArray(gridVAO);
//     glPointSize(5.0f);
//     glDrawArrays(GL_LINES, 0, vertexCount / 3);
//     glBindVertexArray(0);
// }
// std::vector<float> CreateGridVertices(float size, int divisions, const std::vector<Object>& objs) {
    
//     std::vector<float> vertices;
//     float step = size / divisions;
//     float halfSize = size / 2.0f;

//     // x axis
//     for (int yStep = 3; yStep <= 3; ++yStep) {
//         float y = -halfSize*0.3f + yStep * step;
//         for (int zStep = 0; zStep <= divisions; ++zStep) {
//             float z = -halfSize + zStep * step;
//             for (int xStep = 0; xStep < divisions; ++xStep) {
//                 float xStart = -halfSize + xStep * step;
//                 float xEnd = xStart + step;
//                 vertices.push_back(xStart); vertices.push_back(y); vertices.push_back(z);
//                 vertices.push_back(xEnd);   vertices.push_back(y); vertices.push_back(z);
//             }
//         }
//     }
//     // zaxis
//     for (int xStep = 0; xStep <= divisions; ++xStep) {
//         float x = -halfSize + xStep * step;
//         for (int yStep = 3; yStep <= 3; ++yStep) {
//             float y = -halfSize*0.3f + yStep * step;
//             for (int zStep = 0; zStep < divisions; ++zStep) {
//                 float zStart = -halfSize + zStep * step;
//                 float zEnd = zStart + step;
//                 vertices.push_back(x); vertices.push_back(y); vertices.push_back(zStart);
//                 vertices.push_back(x); vertices.push_back(y); vertices.push_back(zEnd);
//             }
//         }
//     }

    

//     return vertices;

// }
// std::vector<float> UpdateGridVertices(std::vector<float> vertices, const std::vector<Object>& objs, float halfSize, float originalY) {

//     glm::vec3 cornerLL(-halfSize, originalY, -halfSize); 
//     float dy_LL = 0.0f;
//     for (const auto& obj : objs) {
//         glm::vec3 toObject = obj.GetPos() - cornerLL;
//         float distance = glm::length(toObject);
//         float distance_m = distance * 1000.0f;
//         float rs = (2 * G * obj.mass) / (c * c);
//         if (distance_m > rs) {
//             float dz = 2 * std::sqrt(rs * (distance_m - rs));
//             dy_LL += dz * 2.0f;
//         }
//     }

//     glm::vec3 cornerLR(halfSize, originalY, -halfSize); 
//     float dy_LR = 0.0f;
//     for (const auto& obj : objs) {
//         glm::vec3 toObject = obj.GetPos() - cornerLR;
//         float distance = glm::length(toObject);
//         float distance_m = distance * 1000.0f;
//         float rs = (2 * G * obj.mass) / (c * c);
//         if (distance_m > rs) {
//             float dz = 2 * std::sqrt(rs * (distance_m - rs));
//             dy_LR += dz * 2.0f;
//         }
//     }

//     glm::vec3 cornerUL(-halfSize, originalY, halfSize); 
//     float dy_UL = 0.0f;
//     for (const auto& obj : objs) {
//         glm::vec3 toObject = obj.GetPos() - cornerUL;
//         float distance = glm::length(toObject);
//         float distance_m = distance * 1000.0f;
//         float rs = (2 * G * obj.mass) / (c * c);
//         if (distance_m > rs) {
//             float dz = 2 * std::sqrt(rs * (distance_m - rs));
//             dy_UL += dz * 2.0f;
//         }
//     }

//     glm::vec3 cornerUR(halfSize, originalY, halfSize); 
//     float dy_UR = 0.0f;
//     for (const auto& obj : objs) {
//         glm::vec3 toObject = obj.GetPos() - cornerUR;
//         float distance = glm::length(toObject);
//         float distance_m = distance * 1000.0f;
//         float rs = (2 * G * obj.mass) / (c * c);
//         if (distance_m > rs) {
//             float dz = 2 * std::sqrt(rs * (distance_m - rs));
//             dy_UR += dz * 2.0f;
//         }
//     }

//     for (size_t i = 0; i < vertices.size(); i += 3) {
//         float x = vertices[i];
//         float z = vertices[i + 2];

//         glm::vec3 vertexPos(x, originalY, z);
//         float dy = 0.0f;
//         for (const auto& obj : objs) {
//             glm::vec3 toObject = obj.GetPos() - vertexPos;
//             float distance = glm::length(toObject);
//             float distance_m = distance * 1000.0f;
//             float rs = (2 * G * obj.mass) / (c * c);
//             if (distance_m > rs) {
//                 float dz = 2 * std::sqrt(rs * (distance_m - rs));
//                 dy += dz * 2.0f;
//             }
//         }

//         float u = (x + halfSize) / (2 * halfSize); 
//         float v = (z + halfSize) / (2 * halfSize);

//         float shift = (1 - u) * (1 - v) * dy_LL +  // Lower-left contribution
//                       u * (1 - v) * dy_LR +        // Lower-right
//                       (1 - u) * v * dy_UL +        // Upper-left
//                       u * v * dy_UR;               // Upper-right

//         vertices[i + 1] = originalY + (dy - shift) + halfSize / 3;
//     }

//     return vertices;
// }













