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
const float electron_r = 0.25f; // Will be used as radius for ray-traced spheres
const double hbar = 1;
const double m_e = 1;
const double zmSpeed = 10.0;

// --- orbital stats ---
int N = 100000;
float LightingScaler = 700;
float n = 3; float l = 1; float m = 1;


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

vec4 inferno2(double r, double theta, double phi, int n, int l, int m)
{
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

    // log compression
    double t = log10(intensity + 1e-12) + 12.0;
    t /= 12.0;

    t = clamp(t, 0.0, 1.0);

    // --- inferno-style ramp ---
    float rC = smoothstep(0.15f, 1.0f, static_cast<float>(t));
    float gC = smoothstep(0.45f, 1.0f, static_cast<float>(t));
    float bC = smoothstep(0.85f, 1.0f, static_cast<float>(t)) * 0.2f;

    return vec4(rC, gC * 0.8f, bC, 1.0f);
}

vec4 heatmap_fire(float value) {
    // Ensure value is clamped between 0 and 1
    value = std::max(0.0f, std::min(1.0f, value));

    // Define color stops for the "Heat/Fire" pattern
    // Order: Black -> Dark Purple -> Red -> Orange -> Yellow -> White
    const int num_stops = 6;
    vec4 colors[num_stops] = {
        {0.0f, 0.0f, 0.0f, 1.0f}, // 0.0: Black
        {0.3f, 0.0f, 0.6f, 1.0f}, // 0.2: Dark Purple
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
    // result = vec4(0.2, 0.9, 0.05, 1.0);

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
    // return vec4(1.0f);
    return heatmap_fire(intensity * LightingScaler); // Scale for better color mapping
}

// ================= Raytracer ================= //
struct Sphere { vec4 center_radius;  vec4 color; };

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

    // Raytracing vals
    GLuint raytracingShaderProgram;
    GLuint fullscreen_VAO, fullscreen_VBO;
    GLuint ssbo_spheres;

    Engine () {
        // --- Init GLFW ---
        if (!glfwInit()) { cerr << "GLFW init failed\n"; exit(EXIT_FAILURE); } 

        // --- Init Window ---
        window = glfwCreateWindow(WIDTH, HEIGHT, "Quantum Simulation by kavan G - Raytraced", nullptr, nullptr);
        if (!window) { cerr << "Failed to create GLFW window\n"; glfwTerminate(); exit(EXIT_FAILURE); } 
        glfwMakeContextCurrent(window); glViewport(0, 0, WIDTH, HEIGHT); glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        
        // --- Init GLEW ---
        glewExperimental = GL_TRUE;
        if (glewInit() != GLEW_OK) { cerr << "Failed to initialize GLEW\n"; glfwTerminate(); exit(EXIT_FAILURE); }

        // blending for smooth rendering
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        // --- NEW: Initialize the SSBO ---
        glGenBuffers(1, &ssbo_spheres);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo_spheres);
        // Allocate space for 25,000 spheres (matches your particle count)
        glBufferData(GL_SHADER_STORAGE_BUFFER, N * sizeof(Sphere), NULL, GL_DYNAMIC_DRAW);
        // Bind the buffer to index 0 so the shader can see it
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, ssbo_spheres);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

        raytracingShaderProgram = CreateRaytracingShaderProgram();
        setupFullscreenQuad();
    }
    vec3 sphericalToCartesian(float r, float theta, float phi){
        float x = r * sin(theta) * cos(phi);
        float y = r * cos(theta);
        float z = r * sin(theta) * sin(phi);
        return vec3(x, y, z);
    }
    
    // ==== Raytracer functions ==== //
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
            uniform float light_intensity;
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
                    vec3 diffuse = diff * sphere_color * shadow_factor * light_intensity;
                    
                    // Ambient lighting
                    vec3 ambient = ambient_light * sphere_color * light_intensity;
                    
                    FragColor = vec4(ambient + diffuse, spheres[closest_sphere_idx].color.a);
                } else {
                    // Nothing hit, draw background
                    FragColor = vec4(0.0f, 0.0f, 0.0f, 0.0f); // Black Background
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
    void runRayTracer(vector<Sphere> sphere_data) {
        if (sphere_data.empty()) return;
        // Update GPU buffer with current sphere positions
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo_spheres);
        glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, sphere_data.size() * sizeof(Sphere), sphere_data.data());
        
        // Crucial: Link the buffer to binding point 0 for the shader
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, ssbo_spheres);

        // Raytracing Render
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glUseProgram(raytracingShaderProgram);

        // Update GPU buffer
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo_spheres);
        glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, sphere_data.size() * sizeof(Sphere), sphere_data.data());
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

        // Raytracing Render
        glClear(GL_COLOR_BUFFER_BIT);
        glUseProgram( raytracingShaderProgram);

        mat4 view = lookAt(camera.position(), camera.target, vec3(0,1,0));
        mat4 projection = perspective(radians(45.0f), (float)WIDTH/HEIGHT, 0.1f, 10000.0f);
        mat4 invViewProj = inverse(projection * view);

        // Define light properties and pass them to the shader
        vec3 light_pos = vec3(0.0f, 50.0f, 50.0f);
        vec3 ambient_light = vec3(0.2f);
        float light_intensity = 3.0f;

        glUniform3fv(glGetUniformLocation( raytracingShaderProgram, "camera_pos"), 1, value_ptr(camera.position()));
        glUniformMatrix4fv(glGetUniformLocation( raytracingShaderProgram, "inv_view_proj"), 1, GL_FALSE, value_ptr(invViewProj));
        glUniform3fv(glGetUniformLocation( raytracingShaderProgram, "light_pos"), 1, value_ptr(light_pos));
        glUniform3fv(glGetUniformLocation( raytracingShaderProgram, "ambient_light"), 1, value_ptr(ambient_light));
        glUniform1f(glGetUniformLocation(raytracingShaderProgram, "light_intensity"), light_intensity);
        
        glBindVertexArray(fullscreen_VAO);
        glDrawArrays(GL_TRIANGLES, 0, 6);
        glBindVertexArray(0);
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
                N *=10;
                generateParticles(N);
            } else if (key == GLFW_KEY_G) {
                N /=10;
                generateParticles(N);
            }

            // Clamp to valid ranges
            if (l > n - 1) l = n - 1;
            if (l < 0) l = 0;
            if (m > l) m = l;
            if (m < -l) m = -l;

            //electron_r = float(n) / 3.0f;
            cout << "Quantum numbers updated: n=" << n << " l=" << l << " m=" << m << " N=" << N << "\n";
        });
    }
};
Engine engine;


// ================= Main Loop ================= //
int main () {
    engine.setupCameraCallbacks();

    // --- Quantum numbers setup ---

    // --- Sample particles ---
    for (int i = 0; i < N; ++i) {
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

    vector<Sphere> spheres;
    for(const auto& p : particles) {
        spheres.push_back({vec4(p.pos, electron_r), p.color});
    }

    float dt = 0.5f;
    cout << "Starting simulation..." << endl;
    while (!glfwWindowShouldClose(engine.window)) {

        // ------ Draw Particles ------
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
            if (p.pos.z < 0 || p.pos.y < 0)
            sphere_data.push_back({vec4(p.pos, electron_r), p.color});
        }
        engine.runRayTracer(sphere_data);

        glfwSwapBuffers(engine.window);
        glfwPollEvents();
    }

    // --- Cleanup ---
    glDeleteBuffers(1, &engine.ssbo_spheres);
    
    glfwDestroyWindow(engine.window);
    glfwTerminate();
    return 0;
}