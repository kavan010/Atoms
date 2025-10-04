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

const float c = 299792458.0f;    // speed of light in m/s
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

// ================= Particles ================= //
struct Proton { float x, y; };
struct Neutron { float x, y; };
struct Electron {float orbitRadius; float angle; float speed; };

struct Atom {
    int protons, neutrons, electrons, x, y;
    vector<Proton> protonList;
    vector<Neutron> neutronList;
    vector<Electron> electronList;

    Atom(int p, int e, int n, int x, int y) : protons(p), neutrons(n), electrons(e), x(x), y(y) {
        srand(time(NULL));

    // Place protons clustered at center
        for (int i = 0; i < protons; i++) {
            float jitterX = (rand() % 10 - 5); // small random offset
            float jitterY = (rand() % 10 - 5);
            protonList.push_back({jitterX + x, jitterY + y});
        }

        // Place neutrons clustered at center
        for (int i = 0; i < neutrons; i++) {
            float jitterX = (rand() % 10 - 5);
            float jitterY = (rand() % 10 - 5);
            neutronList.push_back({jitterX + x, jitterY + y});
        }

        // Distribute electrons into shells (2, 8, 18, 32…)
        int remaining = electrons;
        int shell = 1;
        int shellCap[] = {2, 8, 18, 32}; // enough for simple atoms
        float baseRadius = 45.0f;

        for (int s = 0; remaining > 0 && s < 4; s++) {
            int count = std::min(remaining, shellCap[s]);
            for (int i = 0; i < count; i++) {
                float angle = i * (2.0f * M_PI / count);
                electronList.push_back({baseRadius * (s+1), angle, c / 50000000.0f}); // speed decreases with shell number
            }
            remaining -= count;
        }
    };

    void update() {
        double dt = 1.0 / 60.0; // assuming 60 FPS for simplicity
        for (auto &e : electronList) {
            e.angle += e.speed * dt;
            if (e.angle > 2*M_PI) e.angle -= 2*M_PI;
        }
    }
    void draw() {
        // Draw nucleus (protons red, neutrons gray)
        for (auto &p : protonList) {
            glColor3f(1.0f, 0.0f, 0.0f);
            engine.drawFilledCircle(p.x, p.y, 12);
        }

        for (auto &n : neutronList) {
            glColor3f(0.5f, 0.5f, 0.5f);
            engine.drawFilledCircle(n.x, n.y, 12);
        }

        // Draw electrons (blue)
        glColor3f(0.0f, 0.6f, 1.0f);
        for (auto &e : electronList) {
            float ex = cosf(e.angle) * e.orbitRadius;
            float ey = sinf(e.angle) * e.orbitRadius;
            engine.drawFilledCircle(ex + x, ey + y, 6);
        }
    }
};


vector<Atom> atoms = { 
    //   p  e  n   x    y
    Atom(1, 1, 2, -200, 0),   // Helium
    Atom(6, 6, 12, 200, 0)    // Carbon
};
// ================= Main ================= //
int main () {

    while (!glfwWindowShouldClose(engine.window)) {
        engine.run();

        for (auto &atom : atoms) {
            atom.draw();
            atom.update();
        };

        glfwSwapBuffers(engine.window);
        glfwPollEvents();
    }

    return 0;
}


