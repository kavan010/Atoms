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
float orbitDistance = 15.0f;

// --- engine ---
struct Wave;
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


// --- waves ---
struct WavePoint { vec2 localPos; vec2 dir;  };
struct Wave {
   float energy;
   float sigma = 40.0f, k = 0.4f, phase = 0.0f, a = 10.0f, angleR;
   vector<WavePoint> points;


   vec2 pos, dir;
   vec3 col;


   Wave(float e, vec2 pos, vec2 dir, vec3 col = vec3(0.0f, 1.0f, 1.0f)) : energy(e), pos(pos), dir(dir), col(col) {
       dir = normalize(dir);
       for (float x = -sigma; x <= sigma; x += 0.1f) {
           points.push_back({ pos + x*dir, dir * 200.0f});
       }
       angleR = atan2(dir.y, dir.x);
   }
  
   void draw() {
        glColor3f(col.r, col.g, col.b);
        glBegin(GL_LINE_STRIP);
        for (WavePoint& p : points) {
            // Perpendicular vector for sine displacement
            vec2 perp(-p.dir.y, p.dir.x);
            perp = normalize(perp);


        // Use global phase, not x, for sine
        float y_disp = a * sin(k * length(p.localPos) - phase);


            vec2 drawPos = p.localPos + perp * y_disp;
            glVertex2f(drawPos.x, drawPos.y);
        }
        glEnd();
    }
   bool update(float dt) {
       phase += 30.0f * dt; // continuous phase


       for (WavePoint& p : points) {
           // move along velocity
           p.localPos += p.dir * dt;


            if (p.localPos.x < -engine.WIDTH/2.0f || p.localPos.x > engine.WIDTH/2.0f || p.localPos.y < -engine.HEIGHT/2.0f || p.localPos.y > engine.HEIGHT/2.0f) {
                return true;
            }
       }
       return false;
   }
};
vector<Wave> waves { };


// --- particles ---
struct Particle {
    vec2 pos;
    int charge;
    float angle = 0.0f;
    int n = 1;
    float excitedTimer = 0.0f;
    Particle(vec2 pos, int charge) : pos(pos), charge(charge) {}

    void draw (vec2 centre, int segments = 50) {

        // --- draw outline ---
        if (charge == -1) {
            glLineWidth(0.4f);
            glBegin(GL_LINE_LOOP);
            glColor3f(0.4f, 0.4f, 0.4f);
            for (int i = 0; i <= segments; i++) {
                float angle = 2.0f * M_PI * i/segments;
                float x = cos(angle) * n * orbitDistance;
                float y = sin(angle) * n * orbitDistance;  
                glVertex2f(x + centre.x, y + centre.y);
            }
            glEnd();
        }

        // --- draw particles ---
        float r;
        if (charge == -1)       { r = 2; glColor3f(0.0f, 1.0f, 1.0f); } 
        else if (charge == 1)   { r = 5; glColor3f(1.0f, 0.0f, 0.0f); } 
        else                    { r = 5; glColor3f(0.5f, 0.5f, 0.5f); }

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
        // --- set radius ---
        float r = n * orbitDistance;
        angle += 0.05;
        // --- update pos with new angle and radius ---
        pos = vec2( cos(angle) * r + c.x, 
                    sin(angle) * r + c.y
                );

        if (excitedTimer <= 0.0f && n > 1) {
            n--;
            excitedTimer += 0.003f;
            float waveDirX = ((float)rand() / RAND_MAX) * 2.0f - 1.0f;
            float waveDirY = ((float)rand() / RAND_MAX) * 2.0f - 1.0f;
            float energyDiff = -13.6f/((n+1)*(n+1)) - (-13.6f/(n*n));
            waves.emplace_back(energyDiff, pos, vec2(waveDirX, waveDirY), vec3(1.0f, 1.0f, 0.0f));
        }
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
};



static void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods) {
    if (button != GLFW_MOUSE_BUTTON_LEFT || action != GLFW_PRESS) return;

    double mx, my;
    glfwGetCursorPos(window, &mx, &my);

    Engine* engine = static_cast<Engine*>(glfwGetWindowUserPointer(window));

    // screen → world (matches your glOrtho setup)
    float worldX = (float)mx - engine->WIDTH / 2.0f;
    float worldY = engine->HEIGHT / 2.0f - (float)my;
    vec2 spawnPos(worldX, worldY);

    // spawn 25 waves in all directions
    float energyN1toN2 = -13.6f/(2*2) - (-13.6f);
    for (int i = 0; i < 25; i++) {
        float angle = ((float)rand() / RAND_MAX) * 2.0f * M_PI;
        vec2 dir(cos(angle), sin(angle));

        waves.push_back(
            Wave(energyN1toN2, spawnPos, dir)
        );
    }
}



// --- main ---
int main () {

    // Initialize 20 atoms in a circle at the center
    {
        int num_atoms = 20;
        float radius = 100.0f; // Radius of the circle
        for (int i = 0; i < num_atoms; i++) {
            float angle = 2.0f * M_PI * i / num_atoms;
            float x = cos(angle) * radius;
            float y = sin(angle) * radius;
            atoms.emplace_back(vec2(x, y));
        }
    }
    
    // callbacks
    glfwSetWindowUserPointer(engine.window, &engine);
    glfwSetMouseButtonCallback(engine.window, mouseButtonCallback);

    // -- init waves --
    float energyN1toN2 = -13.6f/(2*2) - (-13.6f);
    for (int i = 0; i < 24; i++) {
        waves.push_back(
            Wave(energyN1toN2, vec2(200.0f, i*20-200), vec2(-1.0f, 0.0f))
        );
    }



    while (!glfwWindowShouldClose(engine.window)) {
        engine.run();

        // --- Draw Particles ---
        
        for (Atom &a : atoms) {
            for (Atom &a2 : atoms) {
                if (&a2 == &a) continue;
                float dist = length(a.pos - a2.pos);
                vec2 dir = normalize(a.pos - a2.pos);
                a.v += dir / dist * 57.5f; // Repulsion force inversely proportional to distance
            }

            // --- Boundary Repulsion ---
            const float boundary_stiffness = 0.01f;
            const float boundary_threshold = 200.0f;

            // Left boundary
            float dist_left = a.pos.x + engine.WIDTH / 2.0f;
            if (dist_left < boundary_threshold) {
                a.v.x += (boundary_threshold - dist_left) * boundary_stiffness;
            }

            // Right boundary
            float dist_right = engine.WIDTH / 2.0f - a.pos.x;
            if (dist_right < boundary_threshold) {
                a.v.x -= (boundary_threshold - dist_right) * boundary_stiffness;
            }

            // Top boundary
            float dist_top = engine.HEIGHT / 2.0f - a.pos.y;
            if (dist_top < boundary_threshold) {
                a.v.y -= (boundary_threshold - dist_top) * boundary_stiffness;
            }

            // Bottom boundary
            float dist_bottom = a.pos.y + engine.HEIGHT / 2.0f;
            if (dist_bottom < boundary_threshold) {
                a.v.y += (boundary_threshold - dist_bottom) * boundary_stiffness;
            }

            
            //a.pos += a.v;
            a.v *= 0.99f; // Damping to stabilize the simulation


            for (Particle &p : a.particles) {
                p.draw(a.pos);

                // --- electrons ---
                if (p.charge == 1) {p.pos = a.pos;}
                if (p.charge == -1) {
                    if (p.excitedTimer > 0.0f) {
                        p.excitedTimer -= 0.001f;
                    }
                    p.update(a.pos);

                    for (Wave& wave : waves) {
                        for (WavePoint& wp : wave.points) {
                            float dist = length(p.pos - wp.localPos);
                            float energyforUp = -13.6f/((p.n+1)*(p.n+1)) - (-13.6f/(p.n*p.n));
                            if (dist < 20.0f && wave.energy == energyforUp && wave.col != vec3(1.0f, 1.0f, 0.0f)) {
                                wave.energy = 0.0f;
                                p.n += 1;
                                p.excitedTimer += 0.003f;
                                break;
                            }
                        }
                    }
                }
            }
        }
               // --- Draw Waves ---
        for (auto it = waves.begin(); it != waves.end(); ) {
            if (it->energy == 0.0f) {
                ++it;
                continue;
            }
            it->draw();
            if (it->update(0.01f)) {
                it = waves.erase(it);
            } else {
                ++it;
            }
        }


        glfwSwapBuffers(engine.window);
        glfwPollEvents();
    }
} 