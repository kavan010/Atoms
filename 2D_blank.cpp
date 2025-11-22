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
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
using namespace glm;
using namespace std;

const float c = 299792458.0f / 100000000.0f;    // speed of light in m/s
const float k = 8.9875517923e9f; // Coulomb's constant

// ================= Engine ================= //
struct Engine {
    GLFWwindow* window;
    int WIDTH = 800;
    int HEIGHT = 600;

    Engine () {
        if(!glfwInit()) {
            cerr << "failed to init glfw, PANIC!";
            exit(EXIT_FAILURE);
        };
        window = glfwCreateWindow(WIDTH, HEIGHT, "Atom Sim", nullptr, nullptr);
        if(!window) {
            cerr << "failed to create window, PANIC!";
            glfwTerminate();
            exit(EXIT_FAILURE);
        }
        glfwMakeContextCurrent(window);
        glViewport(0, 0, WIDTH, HEIGHT);
    };
    void run() {
        glClear(GL_COLOR_BUFFER_BIT);
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        double left   = -WIDTH;
        double right  = WIDTH;
        double bottom = -HEIGHT;
        double top    = HEIGHT;
        glOrtho(left, right, bottom, top, -1.0, 1.0);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
    }
    
    void drawCircle(float x, float y, float r, int segments = 100) {
        glBegin(GL_LINE_LOOP);
        for (int i = 0; i < segments; i++) {
            float angle = 2.0f * M_PI * i / segments;
            float dx = r * cosf(angle);
            float dy = r * sinf(angle);
            glVertex2f(x + dx, y + dy);
        }
        glEnd();
    }
    void drawFilledCircle(float x, float y, float r, int segments = 50) {
        glBegin(GL_TRIANGLE_FAN);
        glVertex2f(x, y); // center
        for (int i = 0; i <= segments; i++) {
            float angle = 2.0f * M_PI * i / segments;
            float dx = r * cosf(angle);
            float dy = r * sinf(angle);
            glVertex2f(x + dx, y + dy);
        }
        glEnd();
    }


};
Engine engine;

// ================= Objects ================= //


// ================= Main ================= //
int main () {

    while (!glfwWindowShouldClose(engine.window)) {
        engine.run();
        

        glfwSwapBuffers(engine.window);
        glfwPollEvents();
    }

    return 0;
}


