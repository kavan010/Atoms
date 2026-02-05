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

// --- variables --- 
float orbitDistance = 200.0f;

// --- engine ---
vec2 mouseWorld(0.0f);
struct Engine {

    GLFWwindow* window;
    int WIDTH = 800, HEIGHT = 600;

    Engine () {
        // --- Init GLFW ---
        if (!glfwInit()) {
            cerr << "failed to init glfw, LOL";
            exit(EXIT_FAILURE);
        }

        // --- Create Window ---
        window = glfwCreateWindow(WIDTH, HEIGHT, "2D atom sim by kavan", nullptr, nullptr);
        if (!window) {
            cerr << "failed to create window, LOLOLOL";
            glfwTerminate();
            exit(EXIT_FAILURE);
        }

        glfwMakeContextCurrent(window);
        int fbWidth, fbHeight;
        glfwGetFramebufferSize(window, &fbWidth, &fbHeight);
        glViewport(0, 0, fbWidth, fbHeight);
    }
    void run() {
        glClear(GL_COLOR_BUFFER_BIT);
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();

        // set origin to centre
        double halfWidth = WIDTH / 2.0f, halfHeight = HEIGHT / 2.0f;
        glOrtho(-halfWidth, halfWidth, -halfHeight, halfHeight, -1.0, 1.0);

        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
    }


};
Engine engine;


// --- particles ---
struct Particle {
    vec2 pos;
    int charge;
    float angle = M_PI;
    float energy = -13.6f;
    int n = 1;
    float excitedTimer = 0.0f;
    Particle(vec2 pos, int charge) : pos(pos), charge(charge) {}

    void draw (vec2 centre, int segments = 50) {

        // --- draw outline ---
        if (charge == -1) {
            segments = 5000;
            glLineWidth(0.4f);
            glBegin(GL_LINE_LOOP);
            glColor3f(0.4f, 0.4f, 0.4f);
            
            float numOsolations = -13.6f / energy;
            float baseOrbit = orbitDistance;
            float amplitude = 50.0f;

            for (int i = 0; i <= segments; i++) {
                float loop_angle = 2.0f * M_PI * i / segments;
                float r = baseOrbit + amplitude * sin(numOsolations * loop_angle);
                float x = cos(loop_angle) * r;
                float y = sin(loop_angle) * r;
                glVertex2f(x + centre.x, y + centre.y);
            }
            glEnd();
        }

        // --- draw particles ---
        float r;
        if (charge == -1)       { r = 10; glColor3f(0.0f, 1.0f, 1.0f); } 
        else if (charge == 1)   { r = 50; glColor3f(1.0f, 0.0f, 0.0f); } 
        else                    { r = 10; glColor3f(0.5f, 0.5f, 0.5f); }

        
        glBegin(GL_TRIANGLE_FAN);
        glVertex2f(pos.x, pos.y);
        for (int i = 0; i <= segments; i++) {
            float angle = 2.0f * M_PI * i/segments;
            float x = cos(angle) * r;
            float y = sin(angle) * r;  
            glVertex2f(x + pos.x, y + pos.y);
        }
        glEnd();

    }
    void update (vec2 c) {
        
        // --- set radius with oscillation ---
        float numOsolations = 0;
        if (energy < 0) {
            numOsolations = -13.6f / energy;
        }
        float baseOrbit = orbitDistance;
        float amplitude = 50.0f;
        float r = baseOrbit + amplitude * sin(numOsolations * angle);

        angle += 0.05;
        
        // --- update pos with new angle and radius ---
        pos = vec2( cos(angle) * r + c.x, 
                    sin(angle) * r + c.y
                );
    }
};

// --- atom struct ---
struct Atom {
    vec2 pos;
    vec2 v = vec2(0.0f);
    vector<Particle> particles = { };
    Atom(vec2 p) : pos(p) {
        particles.emplace_back(pos, 1);                                 // proton
        particles.emplace_back(vec2(pos.x - orbitDistance, pos.y), -1); // electron
    }
};
vector<Atom> atoms {
    Atom(vec2(0.0f, 0.0f))
};

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (action == GLFW_PRESS || action == GLFW_REPEAT)
    {
        if (key == GLFW_KEY_W)
        {
            for (Atom &a : atoms) {
                for (Particle &p : a.particles) {
                    p.energy += 0.01f;
                    p.angle = 0;
                    cout << "Particle energy: " << p.energy << endl;
                }
            }
        }
        else if (key == GLFW_KEY_S)
        {
            for (Atom &a : atoms) {
                for (Particle &p : a.particles) {
                    p.energy -= 0.01f;
                    cout << "Particle energy: " << p.energy << endl;
                }
            }
        }
        if (key == GLFW_KEY_E)
        {
            for (Atom &a : atoms) {
                for (Particle &p : a.particles) {
                    p.energy += 0.1f;
                    cout << "Particle energy: " << p.energy << endl;
                }
            }
        }
        else if (key == GLFW_KEY_D)
        {
            for (Atom &a : atoms) {
                for (Particle &p : a.particles) {
                    p.energy -= 0.1f;
                    cout << "Particle energy: " << p.energy << endl;
                }
            }
        }
        if (key == GLFW_KEY_R)
        {
            for (Atom &a : atoms) {
                for (Particle &p : a.particles) {
                    p.energy += 1.0f;
                    cout << "Particle energy: " << p.energy << endl;
                }
            }
        }
        else if (key == GLFW_KEY_F)
        {
            for (Atom &a : atoms) {
                for (Particle &p : a.particles) {
                    p.energy -= 1.0f;
                    cout << "Particle energy: " << p.energy << endl;
                }
            }
        }
    }
}

// --- main ---
int main () {

    glfwSetKeyCallback(engine.window, key_callback);

    while (!glfwWindowShouldClose(engine.window)) {
        engine.run();

        // --- Draw Particles ---
        
        for (Atom &a : atoms) {

            for (Particle &p : a.particles) {
                p.draw(a.pos);

                if (p.charge == -1)
                    p.update(a.pos);
            }
        }


        glfwSwapBuffers(engine.window);
        glfwPollEvents();
    }
} 