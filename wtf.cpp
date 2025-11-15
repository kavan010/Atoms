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
struct Particle;
static GLuint g_fullscreenVAO = 0;
static void DrawFullscreenQuad() {
    if (g_fullscreenVAO == 0) {
        // Fullscreen triangle (positions + texcoords via vertex shader trick)
        float verts[] = {
            // pos    // tex
            -1.0f, -1.0f,  0.0f, 0.0f,
             3.0f, -1.0f,  2.0f, 0.0f,
            -1.0f,  3.0f,  0.0f, 2.0f
        };
        GLuint vbo;
        glGenVertexArrays(1, &g_fullscreenVAO);
        glGenBuffers(1, &vbo);
        glBindVertexArray(g_fullscreenVAO);
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glBufferData(GL_ARRAY_BUFFER, sizeof(verts), verts, GL_STATIC_DRAW);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)(2 * sizeof(float)));
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindVertexArray(0);
        // keep vbo alive (GL frees with context close)
    }
    glBindVertexArray(g_fullscreenVAO);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
}
// small file reader
static std::string ReadFile(const char* path) {
    std::ifstream in(path, std::ios::in | std::ios::binary);
    if (!in) {
        std::cerr << "Failed to open file: " << path << std::endl;
        return {};
    }
    std::ostringstream contents;
    contents << in.rdbuf();
    in.close();
    return contents.str();
}
static GLuint CompileShader(const char* src, GLenum type) {
    GLuint s = glCreateShader(type);
    glShaderSource(s, 1, &src, nullptr);
    glCompileShader(s);
    GLint ok;
    glGetShaderiv(s, GL_COMPILE_STATUS, &ok);
    if (!ok) {
        GLint len; glGetShaderiv(s, GL_INFO_LOG_LENGTH, &len);
        std::string log(len, '\0');
        glGetShaderInfoLog(s, len, nullptr, &log[0]);
        std::cerr << "Shader compile error:\n" << log << std::endl;
    }
    return s;
}
static GLuint LinkProgram(GLuint vert, GLuint frag) {
    GLuint p = glCreateProgram();
    glAttachShader(p, vert);
    glAttachShader(p, frag);
    glLinkProgram(p);
    GLint ok;
    glGetProgramiv(p, GL_LINK_STATUS, &ok);
    if (!ok) {
        GLint len; glGetProgramiv(p, GL_INFO_LOG_LENGTH, &len);
        std::string log(len, '\0');
        glGetProgramInfoLog(p, len, nullptr, &log[0]);
        std::cerr << "Program link error:\n" << log << std::endl;
    }
    glDetachShader(p, vert);
    glDetachShader(p, frag);
    glDeleteShader(vert);
    glDeleteShader(frag);
    return p;
}

struct SplatRenderer {
    GLuint vao = 0;
    GLuint vbo = 0;
    GLuint program = 0;
    size_t capacity = 0;

    // shader file names expected in working dir
    const char* vertPath = "splat.vert";
    const char* fragPath = "splat.frag";

    void init(size_t maxParticles = 200000) {
        capacity = maxParticles;
        glGenVertexArrays(1, &vao);
        glGenBuffers(1, &vbo);

        glBindVertexArray(vao);
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glBufferData(GL_ARRAY_BUFFER, capacity * 3 * sizeof(float), nullptr, GL_STREAM_DRAW);

        // pos attribute (layout location 0 in splat.vert)
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);

        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindVertexArray(0);

        // compile shaders
        std::string vsrc = ReadFile(vertPath);
        std::string fsrc = ReadFile(fragPath);
        GLuint vs = CompileShader(vsrc.c_str(), GL_VERTEX_SHADER);
        GLuint fs = CompileShader(fsrc.c_str(), GL_FRAGMENT_SHADER);
        program = LinkProgram(vs, fs);
    }

    // upload particle positions each frame
    void updatePositions(const vector<Particle>& atoms) {
        vector<vec3> positions;
        for (const auto& atom : atoms) {
            positions.push_back(atom.position);
        }
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        // Orphan + buffer subdata for streaming
        glBufferData(GL_ARRAY_BUFFER, capacity * 3 * sizeof(float), nullptr, GL_STREAM_DRAW);
        glBufferSubData(GL_ARRAY_BUFFER, 0, positions.size() * sizeof(glm::vec3), positions.data());
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    // draw into current FBO (assumes correct viewport already set)
    void render(const vector<Particle>& atoms,
                const mat4& view, const mat4& proj,
                float pointScreenSize = 40.0f, float intensity = 1.0f, float sigma = 0.25f) {
        vector<vec3> positions;
        for (const auto& atom : atoms) {
            positions.push_back(atom.position);
        }
        if (positions.empty()) return;
        glUseProgram(program);
        GLint loc;

        loc = glGetUniformLocation(program, "view");
        if (loc != -1) glUniformMatrix4fv(loc, 1, GL_FALSE, &view[0][0]);
        loc = glGetUniformLocation(program, "projection");
        if (loc != -1) glUniformMatrix4fv(loc, 1, GL_FALSE, &proj[0][0]);
        loc = glGetUniformLocation(program, "pointScreenSize");
        if (loc != -1) glUniform1f(loc, pointScreenSize);
        loc = glGetUniformLocation(program, "intensity");
        if (loc != -1) glUniform1f(loc, intensity);
        loc = glGetUniformLocation(program, "sigma");
        if (loc != -1) glUniform1f(loc, sigma);

        // update VBO and draw
        //updatePositions(positions);
        glBindVertexArray(vao);
        glEnable(GL_PROGRAM_POINT_SIZE);
        glDrawArrays(GL_POINTS, 0, (GLsizei)positions.size());
        glBindVertexArray(0);
        glUseProgram(0);
    }
};
struct DensityField {
    GLuint fbo = 0;
    GLuint tex = 0;       // single-channel float texture (GL_R16F)
    GLuint pingTex = 0;   // ping-pong for blur
    GLuint pingFBO = 0;
    int width = 0, height = 0;

    // Shaders
    GLuint splatProgram = 0;
    GLuint blurProgram = 0;
    GLuint compositeProgram = 0;

    void init(int w, int h) {
        width = w; height = h;
        // create tex
        glGenTextures(1, &tex);
        glBindTexture(GL_TEXTURE_2D, tex);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_R16F, width, height, 0, GL_RED, GL_FLOAT, nullptr);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

        // ping
        glGenTextures(1, &pingTex);
        glBindTexture(GL_TEXTURE_2D, pingTex);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_R16F, width, height, 0, GL_RED, GL_FLOAT, nullptr);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

        // FBOs
        glGenFramebuffers(1, &fbo);
        glBindFramebuffer(GL_FRAMEBUFFER, fbo);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tex, 0);
        GLenum db = glCheckFramebufferStatus(GL_FRAMEBUFFER);
        if (db != GL_FRAMEBUFFER_COMPLETE) { cerr << "density fbo incomplete\n"; }

        glGenFramebuffers(1, &pingFBO);
        glBindFramebuffer(GL_FRAMEBUFFER, pingFBO);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, pingTex, 0);

        glBindFramebuffer(GL_FRAMEBUFFER, 0);

        // compile shaders (see below)
        splatProgram = CreateSplatProgram();
        blurProgram  = CreateBlurProgram();
        compositeProgram = CreateCompositeProgram();
    }

    void beginCapture() {
        glBindFramebuffer(GL_FRAMEBUFFER, fbo);
        glViewport(0,0, width, height);
        glClearColor(0,0,0,0);
        glClear(GL_COLOR_BUFFER_BIT);
        // additive blending
        glEnable(GL_BLEND);
        glBlendFunc(GL_ONE, GL_ONE);
    }

    void endCapture() {
        glDisable(GL_BLEND);
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
    }

    // call to blur the tex into itself (ping/pong)
    void blur(int iterations = 2) {
        glUseProgram(blurProgram);
        for (int i = 0; i < iterations; ++i) {
            // horizontal pass -> pingFBO
            glBindFramebuffer(GL_FRAMEBUFFER, pingFBO);
            glViewport(0,0,width,height);
            glActiveTexture(GL_TEXTURE0);
            glBindTexture(GL_TEXTURE_2D, tex);
            glUniform2f(glGetUniformLocation(blurProgram,"uDir"), 1.0f/width, 0.0f);
            DrawFullscreenQuad(); // simple VAO that draws a full-screen triangle/quad
            // vertical pass -> fbo
            glBindFramebuffer(GL_FRAMEBUFFER, fbo);
            glActiveTexture(GL_TEXTURE0);
            glBindTexture(GL_TEXTURE_2D, pingTex);
            glUniform2f(glGetUniformLocation(blurProgram,"uDir"), 0.0f, 1.0f/height);
            DrawFullscreenQuad();
        }
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
    }

    // composite to screen, mapping density -> color
    void compositeToScreen(int screenW, int screenH) {
        glUseProgram(compositeProgram);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, tex);
        glViewport(0,0,screenW,screenH);
        DrawFullscreenQuad();
    }

    // utility placeholders you must implement:
    GLuint CreateSplatProgram();
    GLuint CreateBlurProgram();
    GLuint CreateCompositeProgram();
};

SplatRenderer splat;
DensityField density;

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

    vec3 front = normalize(target - position()); // Direction the camera is looking at
    vec3 up = vec3(0.0f, 1.0f, 0.0f); // Up direction for the camera
    float fov = 45.0f; // Field of view in degrees

    mat4 GetViewMatrix() const {
        return lookAt(position(), target, up);
    }

    mat4 GetProjectionMatrix(float aspect) const {
        return perspective(glm::radians(fov), aspect, 0.1f, 1000.0f);
    }
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

        density.init(WIDTH/2, HEIGHT/2);
        splat.init(200000); // capacity, tune as needed

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

struct Particle {
    GLuint VAO, VBO;
    float radius; // in femtometers
    vec4 color; // RGB values between 0 and 1
    vec3 position; // in picometers
    vector<float> vertices = DrawSphere();
    string orbital = "";

    Particle(float r, vec4 col, vec3 pos, string orbit = "") : radius(r), color(col), position(pos), orbital(orbit) {
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



struct Dumbbell {
    GLuint VAO, VBO;
    vec3 position;
    vec4 color;
    vector<float> vertices;
    char axis;
    float max2px;
    float max2py;
    float max2pz;

    // New constructor signature
    Dumbbell(float max2px, float max2py, float max2pz, char axis) : axis(axis), max2px(max2px), max2py(max2py), max2pz(max2pz)
    {
        vertices = DrawDumbell(); // No more arguments needed for DrawDumbell
        engine.CreateVBOVAO(VAO, VBO, vertices.data(), vertices.size());

        position = vec3(0.0f);
        color = vec4(1.0f, 0.5f, 1.0f, 0.1f);
    }
    
    void Draw(GLint objectColorLoc, GLint modelLoc) {
        glUniform4f(objectColorLoc, color.r, color.g, color.b, color.a);
        glm::mat4 model = glm::mat4(1.0f);
        model = glm::translate(model, position);
        glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model));
        glBindVertexArray(VAO);
        glLineWidth(2.0f); // optional: make lines thicker
        glDrawArrays(GL_TRIANGLES, 0, (GLsizei)(vertices.size() / 3)); // draw the stored line segments
        glBindVertexArray(0);
    }
    
    vector<float> DrawDumbell() {
        vector<float> vertices;
        int sectorCount = 36;
        int stackCount = 18;

        float rx = max2px / 2.0f; // radius x
        float ry = max2py;        // radius y
        float rz = max2pz;        // radius z
        float centerX = rx;       // right sphere
        float centerXMirror = -rx; // left sphere

        // store points for one sphere
        vector<vec3> points((stackCount + 1) * (sectorCount + 1));

        // build first sphere
        // vertices.push_back(rx);
        // vertices.push_back(0.0);
        // vertices.push_back(0.0);
        // vertices.push_back(-rx);
        // vertices.push_back(0.0);
        // vertices.push_back(0.0);
        for (int i = 0; i <= stackCount; ++i) {
            float stackAngle = M_PI / 2 - i * (M_PI / stackCount);
            float xy = cosf(stackAngle);
            float y = ry * sinf(stackAngle);

            for (int j = 0; j <= sectorCount; ++j) {
                float sectorAngle = j * (2 * M_PI / sectorCount);
                float x = rx * xy * cosf(sectorAngle) + centerX;
                float z = rz * xy * sinf(sectorAngle);
                points[i * (sectorCount + 1) + j] = vec3(x, y, z);
            }
        }

        // build triangles for right sphere
        for (int i = 0; i < stackCount; ++i) {
            for (int j = 0; j < sectorCount; ++j) {
                int first = i * (sectorCount + 1) + j;
                int second = first + sectorCount + 1;

                // triangle 1
                vertices.push_back(points[first].x);
                vertices.push_back(points[first].y);
                vertices.push_back(points[first].z);

                vertices.push_back(points[second].x);
                vertices.push_back(points[second].y);
                vertices.push_back(points[second].z);

                vertices.push_back(points[first + 1].x);
                vertices.push_back(points[first + 1].y);
                vertices.push_back(points[first + 1].z);

                // triangle 2
                vertices.push_back(points[first + 1].x);
                vertices.push_back(points[first + 1].y);
                vertices.push_back(points[first + 1].z);

                vertices.push_back(points[second].x);
                vertices.push_back(points[second].y);
                vertices.push_back(points[second].z);

                vertices.push_back(points[second + 1].x);
                vertices.push_back(points[second + 1].y);
                vertices.push_back(points[second + 1].z);
            }
        }

        // now build mirrored sphere by flipping X
        for (int i = 0; i < points.size(); ++i) {
            vec3 p = points[i];
            p.x = -p.x; // mirror along x-axis
            points[i] = p;
        }

        // build triangles for mirrored sphere
        for (int i = 0; i < stackCount; ++i) {
            for (int j = 0; j < sectorCount; ++j) {
                int first = i * (sectorCount + 1) + j;
                int second = first + sectorCount + 1;

                // triangle 1
                vertices.push_back(points[first].x);
                vertices.push_back(points[first].y);
                vertices.push_back(points[first].z);

                vertices.push_back(points[second].x);
                vertices.push_back(points[second].y);
                vertices.push_back(points[second].z);

                vertices.push_back(points[first + 1].x);
                vertices.push_back(points[first + 1].y);
                vertices.push_back(points[first + 1].z);

                // triangle 2
                vertices.push_back(points[first + 1].x);
                vertices.push_back(points[first + 1].y);
                vertices.push_back(points[first + 1].z);

                vertices.push_back(points[second].x);
                vertices.push_back(points[second].y);
                vertices.push_back(points[second].z);

                vertices.push_back(points[second + 1].x);
                vertices.push_back(points[second + 1].y);
                vertices.push_back(points[second + 1].z);
            }
        }

        return vertices;
    }

};

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
    if (true) {
        particles.emplace_back(electron_r, vec4(0.0f, 1.0f, 1.0f, 1.0f), electronPos, "1s"); // cyan-ish color
    }
}
void sample2s() {
    float r = sampleR2s();
    float theta = acos(1.0f - 2.0f * static_cast<float>(rand()) / RAND_MAX); // [0, π]
    float phi   = 2.0f * M_PI * static_cast<float>(rand()) / RAND_MAX;       // [0, 2π]
    
    vec3 electronPos = engine.sphericalToCartesian(r, theta, phi);

    // Use a distinct color for visualization, e.g. yellow for 2s
    if (true) {
    particles.emplace_back(electron_r, vec4(0.0f, 1.0f, 1.0f, 1.0f), electronPos, "2s"); // cyan-ish color
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
    if (true) {
        particles.emplace_back(electron_r, vec4(1.0f, 0.0f, 0.0f, 1.0f), electronPos, "2p_x"); // red for 2p_x
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
        particles.emplace_back(electron_r, vec4(1.0f, 0.0f, 1.0f, 1.0f), electronPos, "2p_y"); // red for 2p_y (keeps existing color)
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
    particles.emplace_back(electron_r, vec4(0.0f, 1.0f, 1.0f, 1.0f), electronPos, "2p_z"); // red for 2p_z
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
    particles.emplace_back(electron_r, vec4(1.0f, 0.0f, 0.0f, 1.0f), electronPos, "3p_z"); // cyan-ish for 3p_z
}
// ================= Main ================= //
int main () {
    setupCameraCallbacks(engine.window);
    GLint modelLoc = glGetUniformLocation(engine.shaderProgram, "model");
    GLint objectColorLoc = glGetUniformLocation(engine.shaderProgram, "objectColor");
    glUseProgram(engine.shaderProgram);

    // ---- GENERATE PARTICLES ---- //
    for (int i = 0; i < 2000; ++i) {
        sample1s();
        sample2s();
        sample2p_x();
        sample2p_y();
        sample2p_z();
        sample3p_z();
        sample1s();
        sample1s();
        sample2s();
        sample2s();
        sample2p_x();
        sample2p_x();
    }

    // -------- MAIN LOOP -------- //
    auto lastSampleTime = std::chrono::steady_clock::now();
    const std::chrono::milliseconds sampleInterval(100); // 0.1s

    
    // 2p
    float max2px = 0.0f;
    float max2py = 0.0f;
    float max2pz = 0.0f;


    float maxRad1s = 0.0f;
    float maxRad2s = 0.0f;


    // ------------------ RENDERING LOOP ------------------
    while (!glfwWindowShouldClose(engine.window)) {
        // maxRad1s = 0.0f;
        // maxRad2s = 0.0f;
        engine.run();

        // ---- DRAW GRID ----
        //grid.Draw(objectColorLoc);

        // ---- SAMPLE NEW PARTICLES PERIODICALLY ----
        // sample1s();
        // sample1s();
        // sample2s();
        // sample2s();
        // sample2p_x();
        // sample2p_x();
        // sample2p_y();
        // sample2p_y();

        // ------------------ DRAW PARTICLES ------------------
        for (auto& p : particles) {
            p.Draw(objectColorLoc, modelLoc);
            glDrawArrays(GL_TRIANGLES, 0, p.vertices.size() / 3);
            if (p.orbital == "1s") {
                float r = length(p.position);
                if (r > maxRad1s) {
                    maxRad1s = r;
                }
            }
            if (p.orbital == "2s") {
                float r = length(p.position);
                if (r > maxRad2s) {
                    maxRad2s = r;
                }
            }
            if (p.orbital == "2p_x") {
                if (p.position.x > max2px) {
                    max2px = p.position.x;
                }
                if (p.position.y > max2py) {
                    max2py = p.position.y;
                }
                if (p.position.z > max2pz) {
                    max2pz = p.position.z;
                }
            }
        
        }
        
        // 1) prepare camera matrices used by splat renderer
glm::mat4 view = camera.GetViewMatrix();
glm::mat4 proj = camera.GetProjectionMatrix(float(engine.WIDTH) / float(engine.HEIGHT));

// 2) capture splats into density buffer
density.beginCapture();
// render particle splats into the density FBO
splat.render(particles, view, proj, /*pointScreenSize*/ 48.0f, /*intensity*/ 1.0f, /*sigma*/ 0.3f);
density.endCapture();

// 3) blur the density (separable)
density.blur(2); // iterations (2 is a good start)

// 4) draw your 3D scene normally to default framebuffer
glBindFramebuffer(GL_FRAMEBUFFER, 0);
glViewport(0, 0, engine.WIDTH, engine.HEIGHT);
// ... render atoms, nucleus, other scene geometry here ...

// 5) composite the density texture over the scene
glEnable(GL_BLEND);
glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); // or additive if you prefer
density.compositeToScreen(engine.WIDTH, engine.HEIGHT);
glDisable(GL_BLEND);

// swap buffers as usual



        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glDepthMask(GL_FALSE); // disable depth writes for transparency

        
        // ------------------ DRAW ORBITALS ------------------
        // Particle orbital = Particle(maxRad1s, vec4(1.0f, 1.0f, 0.0f, 0.1f), vec3(0.0f, 0.0f, 0.0f));
        // orbital.Draw(objectColorLoc, modelLoc);
        // glDrawArrays(GL_TRIANGLES, 0, orbital.vertices.size() / 3);
        // glBindVertexArray(0);

        // Particle orbital2 = Particle(maxRad2s, vec4(0.0f, 1.0f, 1.0f, 0.1f), vec3(0.0f, 0.0f, 0.0f));
        // orbital2.Draw(objectColorLoc, modelLoc);
        // glDrawArrays(GL_TRIANGLES, 1, orbital2.vertices.size() / 3);
        // glBindVertexArray(1);

        // Dumbbell dumb(max2px, max2py, max2pz, 'x');
        // dumb.Draw(objectColorLoc, modelLoc);
        //cout<< "max x rad: " << max2px << "   maxyrad:  "<< max2py << "maxzrad:" << max2pz<<endl;
        

        // particles.erase(particles.begin() + 1, particles.end());

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

    
    glfwDestroyWindow(engine.window);
    glfwTerminate();
    return 0;
}