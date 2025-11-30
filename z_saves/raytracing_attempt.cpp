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
#include <nlohmann/json.hpp>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
using namespace glm;
using namespace std;
using json = nlohmann::json;

const float c = 299792458.0f / 100000000.0f;    // speed of light in m/s
const float eu = 2.71828182845904523536f; // Euler's number
const float k = 8.9875517923e9f; // Coulomb's constant
const float a0 = 52.9f; // Bohr radius in pm
const float electron_r = 1.0f;
const float fieldRes = 25.0f;

// ================= Engine ================= //

struct Camera{
    vec3 target;
    float distance;
    float pitch;
    float yaw;
    vec3 position;
    vec3 up;
    
    // mouse handling
    bool middleMousePressed = false;
    double lastX = 0.0, lastY = 0.0;
    float orbitSpeed = 0.4f;
    float zoomSpeed = 20.0f;

    float fov = 60.0f; 
    double lastMovementTime = 0.0;
    // default: look at (0,0,0), 5 units far.
    Camera(vec3 t = vec3(0.0f, 0.0f, -9.0f), float dist = 5.0f, float yawVal = -90.0f, float pitchVal = 0.0f, float fovVal = 90.0f)
        : target(t), distance(dist), yaw(yawVal), pitch(pitchVal), fov(fovVal) {
        up = vec3(0, 1, 0);
        updatePosition();
    }
    void updatePosition() {
        float radYaw = radians(yaw);
        float radPitch = radians(pitch);
        position.x = target.x + distance * cos(radPitch) * cos(radYaw);
        position.y = target.y + distance * sin(radPitch);
        position.z = target.z + distance * cos(radPitch) * sin(radYaw);
    }
    
    // Member function to handle mouse button events.
    void handleMouseButton(int button, int action, int mods, GLFWwindow* window) {
        if (button == GLFW_MOUSE_BUTTON_MIDDLE) {
            if (action == GLFW_PRESS) {
                middleMousePressed = true;
                glfwGetCursorPos(window, &lastX, &lastY);
                lastMovementTime = glfwGetTime();
            } else if (action == GLFW_RELEASE) {
                middleMousePressed = false;
            }
        }
    }
    void handleCursorPosition(double xpos, double ypos, GLFWwindow* window) {
        if (!middleMousePressed)
            return;

        double deltaX = xpos - lastX;
        double deltaY = ypos - lastY;

        // If shift is held, pan the camera's target.
        if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS ||
            glfwGetKey(window, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS) {
            vec3 forward = normalize(target - position);
            vec3 right = normalize(cross(forward, up));
            vec3 camUp = cross(right, forward);
            float panSpeed = 0.005f * distance;
            target += -right * (float)deltaX * panSpeed + camUp * (float)deltaY * panSpeed;
        }
        // Otherwise, orbit the camera.
        else {
            yaw   += (float)deltaX * orbitSpeed;
            pitch += (float)deltaY * orbitSpeed;
            if (pitch > 89.0f)  pitch = 89.0f;
            if (pitch < -89.0f) pitch = -89.0f;
        }
        updatePosition();
        lastX = xpos;
        lastY = ypos;
        lastMovementTime = glfwGetTime();
    }
    void handleScroll(double xoffset, double yoffset, GLFWwindow* window) {
        // If this is the first input, initialize mouse position
        if (lastX == 0 && lastY == 0) {
            glfwGetCursorPos(window, &lastX, &lastY);
        }

        distance -= (float)yoffset * zoomSpeed;
        if (distance < 1.0f)
            distance = 1.0f;
        
        updatePosition();
        lastMovementTime = glfwGetTime();
    }
    static void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods) {
        Camera* cam = static_cast<Camera*>(glfwGetWindowUserPointer(window));
        cam->handleMouseButton(button, action, mods, window);
    }
    static void cursorPositionCallback(GLFWwindow* window, double xpos, double ypos) {
        Camera* cam = static_cast<Camera*>(glfwGetWindowUserPointer(window));
        cam->handleCursorPosition(xpos, ypos, window);
    }
    static void scrollCallback(GLFWwindow* window, double xoffset, double yoffset) {
        Camera* cam = static_cast<Camera*>(glfwGetWindowUserPointer(window));
        cam->handleScroll(xoffset, yoffset, window);
    }

    void registerCallbacks(GLFWwindow* window) {
        glfwSetWindowUserPointer(window, this);
        glfwSetMouseButtonCallback(window, Camera::mouseButtonCallback);
        glfwSetCursorPosCallback(window, Camera::cursorPositionCallback);
        glfwSetScrollCallback(window, Camera::scrollCallback);
    }
};
struct Ray{
    vec3 direction;
    vec3 origin;
    Ray(vec3 o, vec3 d) : origin(o), direction(normalize(d)){}
};
Camera camera;

struct Engine {

    // opengl vars
    GLFWwindow* window;
    GLuint quadVAO;
    GLuint texture;
    GLuint shaderProgram;
    GLuint screenShaderProgram;

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

        // Make context current BEFORE any GL calls that require a context
        glfwMakeContextCurrent(window);

        // now safe to call GL functions
        glViewport(0, 0, WIDTH, HEIGHT);
        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        glEnable(GL_DEPTH_TEST);

        // Enable alpha blending for transparent objects (use standard blend first)
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

        // Create object shader (for 3D objects) and also create a screen shader for the textured quad
        this->shaderProgram = CreateShaderProgram();           // your existing 3D shader
        this->screenShaderProgram = CreateScreenShaderProgram(); // new shader for screen quad

        // create the fullscreen quad VAO and texture and store them
        std::vector<GLuint> vaoTex = QuadVAO(); // returns {VAO, texture}
        this->quadVAO = vaoTex[0];
        this->texture = vaoTex[1];

        // configure texture default storage to avoid GL errors later (optional)
        glBindTexture(GL_TEXTURE_2D, this->texture);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB8, WIDTH, HEIGHT, 0, GL_RGB, GL_UNSIGNED_BYTE, nullptr);
        glBindTexture(GL_TEXTURE_2D, 0);
    }
    void run() {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glUseProgram(shaderProgram);

        mat4 view = lookAt(camera.position, camera.target, vec3(0,1,0));
        mat4 projection = perspective(radians(45.0f), float(WIDTH)/HEIGHT, 0.1f, 10000.0f);

        glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "view"), 1, GL_FALSE, value_ptr(view));
        glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "projection"), 1, GL_FALSE, value_ptr(projection));

        GLint objectColorLoc = glGetUniformLocation(shaderProgram, "objectColor");

        glfwSwapBuffers(window);
        glfwPollEvents();
    }
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
    GLuint CreateScreenShaderProgram() {
        const char* vs = R"(
        #version 330 core
        layout(location = 0) in vec2 aPos;
        layout(location = 1) in vec2 aTexCoord;
        out vec2 TexCoord;
        void main() {
            TexCoord = aTexCoord;
            gl_Position = vec4(aPos.x, aPos.y, 0.0, 1.0);
        }
        )";

        const char* fs = R"(
        #version 330 core
        in vec2 TexCoord;
        out vec4 FragColor;
        uniform sampler2D screenTexture;
        void main() {
            vec3 col = texture(screenTexture, TexCoord).rgb;
            FragColor = vec4(col, 1.0);
        }
        )";

        GLuint v = glCreateShader(GL_VERTEX_SHADER);
        glShaderSource(v, 1, &vs, nullptr);
        glCompileShader(v);
        // TODO: check compile errors

        GLuint f = glCreateShader(GL_FRAGMENT_SHADER);
        glShaderSource(f, 1, &fs, nullptr);
        glCompileShader(f);
        // TODO: check compile errors

        GLuint prog = glCreateProgram();
        glAttachShader(prog, v);
        glAttachShader(prog, f);
        glLinkProgram(prog);
        glDeleteShader(v);
        glDeleteShader(f);
        return prog;
    }


    vector<GLuint> QuadVAO(){
        float quadVertices[] = {
            // positions   // texCoords
            -1.0f,  1.0f,  0.0f, 1.0f,  // top left
            -1.0f, -1.0f,  0.0f, 0.0f,  // bottom left
            1.0f, -1.0f,  1.0f, 0.0f,  // bottom right

            -1.0f,  1.0f,  0.0f, 1.0f,  // top left
            1.0f, -1.0f,  1.0f, 0.0f,  // bottom right
            1.0f,  1.0f,  1.0f, 1.0f   // top right

        };
        
        GLuint VAO, VBO;
        glGenVertexArrays(1, &VAO);
        glGenBuffers(1, &VBO);

        glBindVertexArray(VAO);
        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertices), quadVertices, GL_STATIC_DRAW);

        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)(2 * sizeof(float)));
        glEnableVertexAttribArray(1);

        GLuint texture;
        glGenTextures(1, &texture);
        glBindTexture(GL_TEXTURE_2D, texture);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        std::vector<GLuint> VAOtexture = {VAO, texture};
        return VAOtexture;
    }
    void renderScene(const std::vector<unsigned char>& pixels, int texWidth, int texHeight) {
        // update texture with ray-tracing results
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, texture);
        // Supply the full texture (update)
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB8, texWidth, texHeight, 0, GL_RGB, GL_UNSIGNED_BYTE, pixels.data());
        // set filtering if needed (already set in QuadVAO)

        // clear screen and draw textured quad
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glUseProgram(screenShaderProgram);
        GLint textureLocation = glGetUniformLocation(screenShaderProgram, "screenTexture");
        glUniform1i(textureLocation, 0); // texture unit 0

        glBindVertexArray(quadVAO);
        glDrawArrays(GL_TRIANGLES, 0, 6);
        glBindVertexArray(0);
    };
    vector<int> OptimizeMovement(double lastMovementTime){
        double currentTime = glfwGetTime();
        bool isMoving = (currentTime - lastMovementTime < 0.2);
        int renderFactor = isMoving ? 4 : 2;
        int rWidth = WIDTH / renderFactor;
        int rHeight = HEIGHT / renderFactor;
        std::vector<int> vec = {rWidth, rHeight};
        return vec;
    }

};
Engine engine;

void setupCameraCallbacks(GLFWwindow* window) {
    glfwSetWindowUserPointer(window, &camera);

    glfwSetMouseButtonCallback(window, [](GLFWwindow* win, int button, int action, int mods) {
        Camera* cam = (Camera*)glfwGetWindowUserPointer(win);
        cam->handleMouseButton(button, action, mods, win);
    });

    glfwSetCursorPosCallback(window, [](GLFWwindow* win, double x, double y) {
        Camera* cam = (Camera*)glfwGetWindowUserPointer(win);
        cam->handleCursorPosition(x, y, win);
    });

    glfwSetScrollCallback(window, [](GLFWwindow* win, double xoffset, double yoffset) {
        Camera* cam = (Camera*)glfwGetWindowUserPointer(win);
        cam->handleScroll(xoffset, yoffset, win);
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

struct Particle{
    vec3 pos;
    float radius;

    Particle(vec3 p, float r = 5.0f) : pos(p), radius(r) {}
    // raytracing
    bool Intersect(Ray &ray, float &t){
        vec3 oc = ray.origin - pos;
        float a = glm::dot(ray.direction, ray.direction); // ray direction scale by t
        float b = 2.0f * glm::dot(oc, ray.direction);     // 
        float c = glm::dot(oc, oc) - radius * radius;     // adjustment by sphere radius
        double discriminant = b*b - 4*a*c;
        if(discriminant < 0){return false;}          // no intersection with sphere

        float intercept = (-b - sqrt(discriminant)) / (2.0f*a);
        if(intercept < 0){
            intercept = (-b + sqrt(discriminant)) / (2.0f*a);
            if(intercept<0){return false;}           // intersection is behind origin
        }
        t = intercept;
        return true;
    };
    vec3 getNormal(vec3 &point) const{
        return normalize(point - pos);
    }
};
struct Scene {
    vector<Particle> particles;
    vec3 lightPos;
    
    // Moved light to a better position for shadows
    Scene(vector<Particle> particles) : lightPos(10.0f, 10.0f, 10.0f), particles(particles) {}

    vec3 trace(Ray &ray){
        float closest = INFINITY;
        const Particle* hitObj = nullptr;

        // 1. Find closest intersection
        for(auto& p : particles){
            float t;
            if(p.Intersect(ray, t)){
                if(t < closest) {
                    closest = t;
                    hitObj = &p;
                }
            }
        };

        // 2. Calculate Lighting
        if(hitObj){
            vec3 hitPoint = ray.origin + ray.direction * closest;
            vec3 normal = hitObj->getNormal(hitPoint);
            
            // Calculate Vector towards light
            vec3 lightDirVec = lightPos - hitPoint;
            float lightDist = length(lightDirVec);
            vec3 lightDir = normalize(lightDirVec);

            // --- A. SHADOWS ---
            Ray shadowRay(hitPoint + normal * 0.01f, lightDir); // Offset to prevent "shadow acne"
            bool inShadow = false;
            
            for(auto& obj : particles) {
                float t;
                if(obj.Intersect(shadowRay, t)) {
                    // key fix: Only cast shadow if the object is CLOSER than the light
                    if (t < lightDist) { 
                        inShadow = true;
                        break;
                    }
                }
            }

            // --- B. PHONG LIGHTING MODEL ---
            
            // 1. Ambient (Base color, always visible)
            vec3 objectColor = vec3(0.1f, 0.6f, 0.9f); // Nice Cyan/Blue Electron color
            float ambientStrength = 0.1f;
            vec3 ambient = ambientStrength * objectColor;

            // If in shadow, we only return ambient light
            if (inShadow) {
                return ambient;
            }

            // 2. Diffuse (Matte brightness based on angle)
            float diff = std::max(glm::dot(normal, lightDir), 0.0f);
            vec3 diffuse = diff * objectColor;

            // 3. Specular (Shiny highlight)
            float specularStrength = 0.5f;
            vec3 viewDir = normalize(ray.origin - hitPoint);
            vec3 reflectDir = reflect(-lightDir, normal); 
            // 32 is the "shininess" factor. Higher = smaller, sharper highlight
            float spec = pow(std::max(dot(viewDir, reflectDir), 0.0f), 32); 
            vec3 specular = specularStrength * spec * vec3(1.0f, 1.0f, 1.0f); // White highlight

            // Combine components
            return ambient + diffuse + specular;
        }

        // Background Color (Dark space blue)
        //return vec3(0.05f, 0.05f, 0.1f);
        return vec3(1.0f); 
    }
};

void drawParticle(Particle particle, GLint modelLoc, GLint objectColorLoc) {
    // Draw each particle
    glUniform4f(objectColorLoc, 1.0f, 1.0f, 1.0f, 1.0f); // Red color for particles
    glm::mat4 model = glm::mat4(1.0f);
    model = translate(model, particle.pos);
    model = scale(model, vec3(2.0f)); // Scale to make it visible
    glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model));

    // Draw a simple point for the particle
    glPointSize(1.0f);
    glBegin(GL_POINTS);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glEnd();
}
vector<Particle> LoadWavefunction(const string& filename) {
    vector<Particle> pts;
    std::ifstream file("orbitals/" + filename);
    if (!file.is_open()) {
        cerr << "Failed to open JSON file: " << filename << endl;
        return pts;
    }

    json j;
    file >> j;

    const float bohr_to_pm = 5.29f;

    for (auto& point : j["points"]) {
        float x = point[0].get<float>() * bohr_to_pm;
        float y = point[1].get<float>() * bohr_to_pm;
        float z = point[2].get<float>() * bohr_to_pm;

        // Particle radius (small for electrons)
        float radius = 1.0f;

        // Color: blue for electron
        //vec4 color = vec4(0.2f, 0.5f, 1.0f, 1.0f);
        vec4 color = vec4(1.0f, 0.0f, 1.0f, 1.0f);

        pts.emplace_back(Particle(vec3(x, y, z), electron_r));
    }

    return pts;
}

// ================= Main ================= //
int main () {
    setupCameraCallbacks(engine.window);
    GLint modelLoc = glGetUniformLocation(engine.shaderProgram, "model");
    GLint objectColorLoc = glGetUniformLocation(engine.shaderProgram, "objectColor");
    glUseProgram(engine.shaderProgram);

    
    vector<Particle> particles = 
    //{
    //     Particle(vec3(0.0f, 0.0f, 0.0f), 10.0f),
    //     Particle(vec3(20.0f, 0.0f, 0.0f), 5.0f),
    //     Particle(vec3(-20.0f, 0.0f, 0.0f), 5.0f)
    // };
    
    LoadWavefunction("orbital_n6_l4_m1.json");
    Scene scene(particles);

    // ------------------ RENDERING LOOP ------------------
    while (!glfwWindowShouldClose(engine.window)) {
        //engine.run();

        int rWidth = 800;
        int rHeight = 600;
        vector<unsigned char> pixels(rWidth * rHeight * 3);
        // render texture (pxl by pxl)
        for(int y = 0; y < rHeight; ++y){
            for(int x = 0; x < rWidth; ++x){
                float scale = tan(radians(camera.fov * 0.5f));

                float aspectRatio = float(rWidth) / float(rHeight);
                float u = float(x) / float(rWidth);
                float v = float(y) / float(rHeight);
                
                // Convert screen coordinates to camera space coordinates with FOV adjustment
                float x_camera = (2.0f * u - 1.0f) * aspectRatio * scale;
                float y_camera = (1.0f - 2.0f * v) * scale; // (1 - 2*v) is equivalent to -(2*v - 1)

                // Transform the ray from camera space to world space.
                vec3 forward = normalize(camera.target - camera.position);
                vec3 right = normalize(cross(forward, vec3(0.0f, 1.0f, 0.0f)));
                vec3 up = cross(right, forward);

                vec3 direction = normalize(x_camera * right + y_camera * up + forward);

                Ray ray(camera.position, direction);
                vec3 color = scene.trace(ray);

                int index = (y * rWidth  + x) * 3;
                pixels[index + 0] = static_cast<unsigned char>(color.r * 255);
                pixels[index + 1] = static_cast<unsigned char>(color.g * 255);
                pixels[index + 2] = static_cast<unsigned char>(color.b * 255);
            }
        }

        // ---- 2. Draw particles -------------------
        // for (const auto& particle : particles) {
        //     drawParticle(particle, modelLoc, objectColorLoc);
        // }
        engine.renderScene(pixels, rWidth, rHeight);
        glfwSwapBuffers(engine.window);
        glfwPollEvents();
    }
    
    glfwDestroyWindow(engine.window);
    glfwTerminate();
    return 0;
}