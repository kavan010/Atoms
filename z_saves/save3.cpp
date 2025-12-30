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
const float a0 = 1;
const float electron_r = 0.3f; // Will be used as radius for ray-traced spheres
const double hbar = 1;
const double m_e = 1;
const double zmSpeed = 10.0;

// ================= Raytracer ================= //
struct Sphere {
    vec4 center_radius; 
    vec4 color;
};

// ================= Engine ================= //
struct Camera {
    vec3 target = vec3(0.0f, 0.0f, 0.0f);
    float radius = 500.0f;
    float azimuth = 0.0f;
    float elevation = M_PI / 2.0f;
    float orbitSpeed = 0.01f;
    float panSpeed = 0.01f;
    double zoomSpeed = zmSpeed;
    bool dragging = false;
    bool panning = false;
    double lastX = 0.0, lastY = 0.0;

    vec3 position() const {
        float clampedElevation = glm::clamp(elevation, 0.01f, float(M_PI) - 0.01f);
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
    
    // Raytracing resources
    GLuint raytracingShaderProgram;
    GLuint fullscreen_VAO, fullscreen_VBO;
    GLuint ssbo_spheres;

    Engine () {
        if (!glfwInit()) { cerr << "GLFW init failed\n"; exit(EXIT_FAILURE); } 

        window = glfwCreateWindow(WIDTH, HEIGHT, "Quantum Simulation by kavan G - Raytraced", nullptr, nullptr);
        if (!window) { cerr << "Failed to create GLFW window\n"; glfwTerminate(); exit(EXIT_FAILURE); } 
        
        glfwMakeContextCurrent(window);
        glViewport(0, 0, WIDTH, HEIGHT);
        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        
        glewExperimental = GL_TRUE;
        if (glewInit() != GLEW_OK) { cerr << "Failed to initialize GLEW\n"; glfwTerminate(); exit(EXIT_FAILURE); }

        // Blending for smooth edges from raytracer
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        raytracingShaderProgram = CreateRaytracingShaderProgram();
        setupFullscreenQuad();
    }

    void setupFullscreenQuad() {
        float quadVertices[] = {
            -1.0f,  1.0f,
            -1.0f, -1.0f,
             1.0f, -1.0f,
            -1.0f,  1.0f,
             1.0f, -1.0f,
             1.0f,  1.0f
        };
        glGenVertexArrays(1, &fullscreen_VAO);
        glGenBuffers(1, &fullscreen_VBO);
        glBindVertexArray(fullscreen_VAO);
        glBindBuffer(GL_ARRAY_BUFFER, fullscreen_VBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertices), &quadVertices, GL_STATIC_DRAW);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)0);
        glBindVertexArray(0);
    }

    GLuint CreateRaytracingShaderProgram() {
        const char* vertexShaderSource = R"(
            #version 430 core
            layout (location = 0) in vec2 aPos;
            out vec2 ScreenPos;
            void main() {
                ScreenPos = aPos;
                gl_Position = vec4(aPos.x, aPos.y, 0.0, 1.0);
            }
        )";

        const char* fragmentShaderSource = R"(
            #version 430 core
            out vec4 FragColor;
            in vec2 ScreenPos;

            uniform vec3 camera_pos;
            uniform mat4 inv_view_proj;
            uniform vec3 light_pos;
            uniform vec3 ambient_light;

            struct Sphere {
                vec4 center_radius; // xyz = center, w = radius
                vec4 color;
            };

            layout(std430, binding = 0) buffer SphereBuffer {
                Sphere spheres[];
            };

            // Ray-sphere intersection test. Returns distance t or -1.0 for no hit.
            float intersect_sphere(vec3 ray_origin, vec3 ray_dir, vec3 sphere_center, float sphere_radius) {
                vec3 oc = ray_origin - sphere_center;
                float a = dot(ray_dir, ray_dir);
                float b = 2.0 * dot(oc, ray_dir);
                float c = dot(oc, oc) - sphere_radius * sphere_radius;
                float discriminant = b * b - 4.0 * a * c;
                if (discriminant < 0.0) {
                    return -1.0;
                } else {
                    return (-b - sqrt(discriminant)) / (2.0 * a);
                }
            }

            // Checks for any intersection along a ray up to a max distance. For shadows.
            bool any_hit(vec3 ray_origin, vec3 ray_dir, float max_dist) {
                uint num_spheres = spheres.length();
                for (uint i = 0; i < num_spheres; ++i) {
                    float t = intersect_sphere(ray_origin, ray_dir, spheres[i].center_radius.xyz, spheres[i].center_radius.w);
                    if (t > 0.0 && t < max_dist) {
                        return true;
                    }
                }
                return false;
            }

            void main() {
                // Construct ray from camera through the screen
                vec4 target = inv_view_proj * vec4(ScreenPos, 1.0, 1.0);
                vec3 ray_dir = normalize(vec3(target / target.w) - camera_pos);

                // Find closest sphere intersection
                float t_min = 1e20;
                int closest_sphere_idx = -1;
                uint num_spheres = spheres.length();
                for (uint i = 0; i < num_spheres; ++i) {
                    float t = intersect_sphere(camera_pos, ray_dir, spheres[i].center_radius.xyz, spheres[i].center_radius.w);
                    if (t > 0.0 && t < t_min) {
                        t_min = t;
                        closest_sphere_idx = int(i);
                    }
                }

                if (closest_sphere_idx != -1) {
                    // We hit a sphere, calculate lighting
                    vec3 hit_pos = camera_pos + t_min * ray_dir;
                    vec3 normal = normalize(hit_pos - spheres[closest_sphere_idx].center_radius.xyz);
                    vec3 sphere_color = spheres[closest_sphere_idx].color.rgb;

                    // Shadow calculation
                    vec3 light_dir = normalize(light_pos - hit_pos);
                    float light_dist = length(light_pos - hit_pos);
                    float shadow_factor = 1.0;
                    // Cast shadow ray, offsetting origin slightly to avoid self-shadowing
                    if (any_hit(hit_pos + normal * 0.001, light_dir, light_dist)) {
                        shadow_factor = 0.0; // Point is in shadow
                    }

                    // Diffuse lighting from point light
                    float diff = max(dot(normal, light_dir), 0.0);
                    vec3 diffuse = diff * sphere_color * shadow_factor;
                    
                    // Ambient lighting
                    vec3 ambient = ambient_light * sphere_color;
                    
                    FragColor = vec4(ambient + diffuse, spheres[closest_sphere_idx].color.a);
                } else {
                    // Nothing hit, draw background
                    FragColor = vec4(1.0, 1.0, 1.0, 1.0); // White Background
                }
            }
        )";

        GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
        glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
        glCompileShader(vertexShader);
        // Error checking here...

        GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
        glShaderSource(fragmentShader, 1, &fragmentShaderSource, NULL);
        glCompileShader(fragmentShader);
        // Error checking here...

        GLuint shaderProgram = glCreateProgram();
        glAttachShader(shaderProgram, vertexShader);
        glAttachShader(shaderProgram, fragmentShader);
        glLinkProgram(shaderProgram);
        // Error checking here...

        glDeleteShader(vertexShader);
        glDeleteShader(fragmentShader);

        return shaderProgram;
    }

    vec3 sphericalToCartesian(float r, float theta, float phi){
        float x = r * sin(theta) * cos(phi);
        float y = r * cos(theta);
        float z = r * sin(theta) * sin(phi);
        return vec3(x, y, z);
    }
};
Engine engine;

void setupCameraCallbacks(GLFWwindow* window) {
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

// ================= Objects ================= //
struct Particle {
    vec3 pos;
    vec3 vel = vec3(0.0f);
    vec4 color = vec4(0.2f, 0.8f, 1.0f, 0.8f); // Corrected alpha
    Particle(vec3 p) : pos(p){}
};
vector<Particle> particles;

// ================= Physics & Sampling ================= //
random_device rd;
mt19937 gen(rd());
uniform_real_distribution<float> dis(0.0f, 1.0f);

double sampleR(int n, int l, mt19937& gen) {
    const int N = 4096;
    const double rMax = 10.0 * n * n * a0;
    static vector<double> cdf;
    if (cdf.empty()) {
        cdf.resize(N);
        double dr = rMax / (N - 1);
        double sum = 0.0;
        for (int i = 0; i < N; ++i) {
            double r = i * dr;
            double rho = 2.0 * r / (n * a0);
            double L = std::assoc_laguerre(n - l - 1, 2 * l + 1, rho);
            double R_norm = sqrt(pow(2.0 / (n * a0), 3) * tgamma(n - l) / (2.0 * n * tgamma(n + l + 1)));
            double R = R_norm * exp(-rho / 2.0) * pow(rho, l) * L;
            double pdf = r * r * R * R;
            sum += pdf;
            cdf[i] = sum;
        }
        for (double& v : cdf) v /= sum;
    }
    uniform_real_distribution<double> u_dis(0.0, 1.0);
    double u = u_dis(gen);
    int idx = lower_bound(cdf.begin(), cdf.end(), u) - cdf.begin();
    return idx * (rMax / (N - 1));
}
double sampleTheta(int l, int m, mt19937& gen) {
    const int N = 2048;
    static vector<double> cdf;
    if (cdf.empty()) {
        cdf.resize(N);
        double dtheta = M_PI / (N - 1);
        double sum = 0.0;
        for (int i = 0; i < N; ++i) {
            double theta = i * dtheta;
            double x = cos(theta);
            double Plm = std::sph_legendre(l, m, x);
            double pdf = sin(theta) * Plm * Plm;
            sum += pdf;
            cdf[i] = sum;
        }
        for (double& v : cdf) v /= sum;
    }
    uniform_real_distribution<double> u_dis(0.0, 1.0);
    double u = u_dis(gen);
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
    double sinTheta = sin(theta);  if (abs(sinTheta) < 1e-4) sinTheta = 1e-4;
    double v_mag = hbar * m / (m_e * r * sinTheta);
    double vx = -v_mag * sin(phi);
    double vy = 0.0;
    double vz =  v_mag * cos(phi);
    return vec3((float)vx, (float)vy, (float)vz);
}

// ================= Main ================= //
int main () {
    setupCameraCallbacks(engine.window);

    float n = 4; float l = 2; float m = 1;

    cout << "Creating particles..." << endl;
    for (int i = 0; i < 15000; ++i) { // Reduced particles for performance
        float r = sampleR(n, l, gen);
        float theta = sampleTheta(l, m, gen);
        float phi = samplePhi(n, l, m);
        vec3 pos = engine.sphericalToCartesian(r, theta, phi);
        particles.emplace_back(pos);
    }
    cout << "Particles created. Setting up render buffer." << endl;

    vector<Sphere> spheres;
    for(const auto& p : particles) {
        spheres.push_back({vec4(p.pos, electron_r), p.color});
    }

    glGenBuffers(1, &engine.ssbo_spheres);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, engine.ssbo_spheres);
    glBufferData(GL_SHADER_STORAGE_BUFFER, spheres.size() * sizeof(Sphere), spheres.data(), GL_DYNAMIC_DRAW);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, engine.ssbo_spheres);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

    float dt = 0.5f;
    cout << "Starting simulation..." << endl;
    while (!glfwWindowShouldClose(engine.window)) {
        // Physics Update
        vector<Sphere> sphere_data;
        for (Particle& p : particles) {
            double r = length(p.pos);
            if (r > 1e-6) {
                double theta = acos(p.pos.y / r);
                p.vel = calculateProbabilityFlow(p, n, l, m);
                vec3 temp_pos = p.pos + p.vel * dt;
                double new_phi = atan2(temp_pos.z, temp_pos.x);
                p.pos = engine.sphericalToCartesian(r, theta, new_phi);
            }
            sphere_data.push_back({vec4(p.pos, electron_r), p.color});
        }
        
        // Update GPU buffer
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, engine.ssbo_spheres);
        glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, sphere_data.size() * sizeof(Sphere), sphere_data.data());
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

        // Raytracing Render
        glClear(GL_COLOR_BUFFER_BIT);
        glUseProgram(engine.raytracingShaderProgram);

        mat4 view = lookAt(camera.position(), camera.target, vec3(0,1,0));
        mat4 projection = perspective(radians(45.0f), (float)engine.WIDTH/engine.HEIGHT, 0.1f, 10000.0f);
        mat4 invViewProj = inverse(projection * view);

        // Define light properties and pass them to the shader
        vec3 light_pos = vec3(200.0f, 300.0f, 200.0f);
        vec3 ambient_light = vec3(0.2f);

        glUniform3fv(glGetUniformLocation(engine.raytracingShaderProgram, "camera_pos"), 1, value_ptr(camera.position()));
        glUniformMatrix4fv(glGetUniformLocation(engine.raytracingShaderProgram, "inv_view_proj"), 1, GL_FALSE, value_ptr(invViewProj));
        glUniform3fv(glGetUniformLocation(engine.raytracingShaderProgram, "light_pos"), 1, value_ptr(light_pos));
        glUniform3fv(glGetUniformLocation(engine.raytracingShaderProgram, "ambient_light"), 1, value_ptr(ambient_light));
        
        glBindVertexArray(engine.fullscreen_VAO);
        glDrawArrays(GL_TRIANGLES, 0, 6);
        glBindVertexArray(0);

        glfwSwapBuffers(engine.window);
        glfwPollEvents();
    }
    
    glDeleteBuffers(1, &engine.ssbo_spheres);
    glfwDestroyWindow(engine.window);
    glfwTerminate();
    return 0;
}