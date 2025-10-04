#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <vector>
#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <sstream>
#include <iomanip>
#include <cstring>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
using namespace glm;
using namespace std;

double G = 6.67430e-11;
double c = 299792458.0;

struct Ray;
void rk4Step (Ray& ray, double dλ, double r_s);

struct Engine {
    GLFWwindow* window;
    int WIDTH = 800;
    int HEIGHT = 600;
    float width = 1.0e11;
    float height = 7.5e10;

    Engine () {
        if(!glfwInit()) {
            cerr << "failed to init glfw, PANIC!";
            exit(EXIT_FAILURE);
        };
        window = glfwCreateWindow(WIDTH, HEIGHT, "2D Black Hole Simulation", nullptr, nullptr);
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
        double left   = -width;
        double right  = width;
        double bottom = -height;
        double top    = height;
        glOrtho(left, right, bottom, top, -1.0, 1.0);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
    }
};
Engine engine;

struct BlackHole {
    vec2 position;
    double mass;
    double r_s;

    BlackHole (vec2 pos, double m) : position(pos), mass(m) {r_s = (2.0 * G * mass) / (c*c);}

    void draw() {
        glColor3f(1.0f, 0.0f, 0.0f); // color of black hole - Red
        glBegin(GL_TRIANGLE_FAN);
        glVertex2f(position.x, position.y);

        for (int i = 0; i <= 100; i++) {
            float angle = 2.0f * M_PI * i / 100;
            float x = cos(angle) * r_s + position.x;
            float y = sin(angle) * r_s + position.y;
            glVertex2f(x, y);
        }

        glEnd();
    }
};
BlackHole SagA(vec2(0.0f, 0.0f), 8.54e36);

struct Ray {
    // -- cartessian coords -- //
    double x; double y;
    // -- polar coords -- //
    double r;   double phi;
    double dr;  double dphi;
    double d2r; double d2phi;
    double E, L;

    vec2 dir;
    vector<vec2> trail;

    Ray (vec2 pos, vec2 dir) : x(pos.x), y(pos.y), dir(dir) {
        r   = hypot(x, y);
        phi = atan(y, x);
        dr   = dir.x * c * cos(phi) + dir.y * c * sin(phi);
        dphi = ( -dir.x * c * sin(phi) + dir.y * c * cos(phi) ) / r;

        L = r*r * dphi;
        double f = 1.0 - SagA.r_s/r;  
        double dt_dλ = sqrt( (dr*dr)/(f*f) + (r*r*dphi*dphi)/f );
        E = f * dt_dλ;

        trail.push_back({x, y});
    };

    void draw () {
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glLineWidth(1.0f);

        size_t N = trail.size();
        glBegin(GL_LINE_STRIP);
        for (size_t i = 0; i < N; ++i) {
            if (N<2) continue;
            float alpha = float(i) / float(N-1);
            glColor4f(1.0f, 1.0f, 1.0f, std::max(alpha, 0.05f));
            glVertex2f(trail[i].x, trail[i].y);
        }
        glEnd();
    }

    void step (float dλ, double r_s) {
        if (r < r_s) return;

        rk4Step(*this, dλ, r_s);

        x = cos(phi) * r;
        y = sin(phi) * r;
        trail.push_back({float(x), float(y)});
    }
};
vector<Ray> rays;

void Geodesic (Ray &ray, double rhs[4], double r_s) {
    float r    = ray.r;
    float phi  = ray.phi; 
    float dr   = ray.dr;
    float dphi = ray.dphi;
    double E    = ray.E;

    double f = 1.0 - r_s/r;


    rhs[0] = dr;
    rhs[1] = dphi;
    double dt_dλ = E / f;
    rhs[2] = 
        - (r_s/(2*r*r)) * f * (dt_dλ*dt_dλ)
        + (r_s/(2*r*r*f)) * (dr*dr)
        + (r - r_s) * (dphi*dphi);
    rhs[3] = -2.0 * dr * dphi / r;
}
void addState (const double a[4], const double b[4], double factor, double out[4]) {
    for (int i = 0; i < 4; i++) {
        out[i] = a[i] + b[i] * factor;
    }
}
void rk4Step (Ray& ray, double dλ, double r_s) {
    double y0[4] = { ray.r, ray.phi, ray.dr, ray.dphi };
    double k1[4], k2[4], k3[4], k4[4], temp[4];

    Geodesic(ray, k1, r_s);
    addState(y0, k1, dλ/2.0, temp);
    Ray r2 = ray; r2.r=temp[0]; r2.phi=temp[1]; r2.dr=temp[2]; r2.dphi=temp[3];
    Geodesic(r2, k2, r_s);

    addState(y0, k2, dλ/2.0, temp);
    Ray r3 = ray; r3.r=temp[0]; r3.phi=temp[1]; r3.dr=temp[2]; r3.dphi=temp[3];
    Geodesic(r3, k3, r_s);

    addState(y0, k3, dλ/2.0, temp);
    Ray r4 = ray; r4.r=temp[0]; r4.phi=temp[1]; r4.dr=temp[2]; r4.dphi=temp[3];
    Geodesic(r4, k4, r_s);

    ray.r    += (dλ/6.0)*(k1[0] + 2*k2[0] + 2*k2[0] + k4[0]);
    ray.phi  += (dλ/6.0)*(k1[1] + 2*k2[1] + 2*k2[1] + k4[1]);
    ray.dr   += (dλ/6.0)*(k1[2] + 2*k2[2] + 2*k2[2] + k4[2]);
    ray.dphi += (dλ/6.0)*(k1[3] + 2*k2[3] + 2*k2[3] + k4[3]);
}

int main () {
    for (float y = -engine.height; y <= engine.height; y+=1e10) {
        rays.push_back(Ray(vec2(-1e11, y), vec2(1.0f, 0.0f)));
    }

    while (!glfwWindowShouldClose(engine.window)) {
        engine.run();
        SagA.draw();

        for (auto& ray : rays) {
            ray.step(1e0, SagA.r_s);
            ray.draw();
        }

        glfwSwapBuffers(engine.window);
        glfwPollEvents();
    }

    return 0;
}


