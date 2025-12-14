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


float c = 299792458.0f;
float h = 6.626e-34f;   // Planck constant
float orbitalRad = 25.0f;
float e = 2.71828182845904523536f;


// ================= Engine ================= //
struct Engine {
   GLFWwindow* window;
   int WIDTH = 800;
   int HEIGHT = 600;

   Engine () {

       // --- init GLFW ---
       if(!glfwInit()) {
           cerr << "failed to init glfw, PANIC!";
           exit(EXIT_FAILURE);
       };

       // --- create window ---
       window = glfwCreateWindow(WIDTH, HEIGHT, "Atom Sim", nullptr, nullptr);
       if(!window) {
           cerr << "failed to create window, PANIC!";
           glfwTerminate();
           exit(EXIT_FAILURE);
       }
       glfwMakeContextCurrent(window);
       int fbWidth, fbHeight;
       glfwGetFramebufferSize(window, &fbWidth, &fbHeight);
       glViewport(0, 0, fbWidth, fbHeight);


   };
   void run() {
       glClear(GL_COLOR_BUFFER_BIT);
       glMatrixMode(GL_PROJECTION);
       glLoadIdentity();


       // Set origin (0,0) to center
       double halfWidth  = WIDTH / 2.0;
       double halfHeight = HEIGHT / 2.0;


       glOrtho(-halfWidth, halfWidth, -halfHeight, halfHeight, -1.0, 1.0);


       glMatrixMode(GL_MODELVIEW);
       glLoadIdentity();
   }
  
   // ----- Drawing Functions -----
   void drawCircle(float x, float y, float r, int segments = 100) {
       glLineWidth(0.4f);
       glColor3f(0.4f, 0.4f, 0.4f);
       glBegin(GL_LINE_LOOP);
       for (int i = 0; i < segments+1; i++) {
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
   float timeDT(float* lastTime) {
       float currentTime = glfwGetTime();
       float dt = currentTime - *lastTime;
       *lastTime = currentTime;
       return dt;
   }


};
Engine engine;

float random(float min, float max, float step) {
   int numSteps = static_cast<int>((max - min) / step + 0.5f);
   int r = std::rand() % (numSteps + 1); // random step index
   return min + r * step;
}
// ================= Particles ================= //


struct WavePoint { vec2 localPos; vec2 dir;  };
struct Wave {
   float energy;
   float sigma = 40.0f, k = 0.4f, phase = 0.0f, a = 10.0f, angleR;
   vector<WavePoint> points;


   vec2 pos, dir;


   Wave(float e, vec2 pos, vec2 dir) : energy(e), pos(pos), dir(dir) {
       dir = normalize(dir);
       for (float x = -sigma; x <= sigma; x += 0.1f) {
           points.push_back({ pos + x*dir, dir * 200.0f});
       }
       angleR = atan2(dir.y, dir.x);
   }
  
   void Draw() {
    
        if (energy == 4.2f) glColor3f(1.0f, 1.0f, 0.0f);
        else glColor3f(0.0f, 1.0f, 1.0f);
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
   void update(float dt) {
       phase += 30.0f * dt; // continuous phase


       for (WavePoint& p : points) {
           // move along velocity
           p.localPos += p.dir * dt;


           // bounce off walls
           if (p.localPos.x < -400) { p.localPos.x = -400; p.dir.x *= -1; }
           if (p.localPos.x > 400)  { p.localPos.x = 400;  p.dir.x *= -1; }
           if (p.localPos.y < -300) { p.localPos.y = -300; p.dir.y *= -1; }
           if (p.localPos.y > 300)  { p.localPos.y = 300;  p.dir.y *= -1; }
       }
   }
};


struct Particle {
   int x, y, charge;
   float theta;

   int n = 1;                      // quantum level
   int targetN = 1;                // excited state target
   bool excited = false;
   float exciteTimer = 0.0f;       // time spent excited
   float exciteLifetime = 0.0f;    // random lifetime

   Particle (int x, int y, int c) : x(x), y(y), charge(c) { theta = atan2f(y, x); }

   void Draw() {
       glColor4f(1.0f, 1.0f, 1.0f, 0.35f);
       engine.drawCircle(0.0f, 0.0f, length(vec2((float) x, (float) y)), 100);

       if (charge == 1){
           glColor3f(1.0f, 0.0f, 0.0f);
           engine.drawFilledCircle(x, y, 10);
       }
       if (charge == 0) {
           glColor3f(0.5f, 0.5f, 0.5f);
           engine.drawFilledCircle(x, y, 10);
       }
       if (charge == -1) {
           glColor3f(0.0f, 1.0f, 1.0f);
           engine.drawFilledCircle(x, y, 2);
       }
   }
   void electronUpdate(float speed, vector<Wave>& waves) {
       theta += speed;
       x = (orbitalRad * (n*n)) * cosf(theta);
       y = (orbitalRad * (n*n)) * sinf(theta);
       if (n!= 1 && exciteTimer <=0.0f) {
           n-=1;
           exciteTimer = 0.003f;
           float waveDirX = (rand() % 2 == 0) ? -1 : 1;
           float waveDirY = (rand() % 2 == 0) ? -1 : 1;
           waves.emplace_back(4.2f, vec2((float)x, (float)y), vec2(waveDirX, waveDirY));
       }
   }
   void checkAbsorbtion(Wave& wave) {
       for (WavePoint& wp : wave.points) {
            float dist = length(vec2((float)x, (float)y) - wp.localPos);
            if (dist < 20.0f && wave.energy == 3.2f) {
                wave.energy = 0.0f;
                n += 1;
                exciteTimer += 0.003f;
                exciteLifetime = random(1.0f, 3.0f, 0.1f);
                break;
            }
        }

   }
};
// ================= Declare Atom Objects ================= //
vector<Particle> particles = {
   Particle(0, 0, 1),
   Particle(50, 50, -1),
};


vector<Wave> waves = { };
// ================= Main ================= //
int main () {
   srand(static_cast<unsigned int>(std::time(nullptr)));


   for (int i = 0; i < 25; i++) {
       float x = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/800.0f)) - 400.0f;
       float y = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/600.0f)) - 300.0f;
       waves.emplace_back(3.2f, vec2(x, y), vec2(y, x));
   }
  
   float lastTime = glfwGetTime(); float dt;
   while (!glfwWindowShouldClose(engine.window)) {

       // --- Time Calculation ---
       dt = engine.timeDT(&lastTime);

       engine.run();

       // ---- Draw Particles ----
       for (Particle& p : particles) {
           if (p.exciteTimer > 0) p.exciteTimer -= dt/100;

           p.Draw();
           if (p.charge == -1){
               p.electronUpdate(5.0f * dt, waves);
               cout<< "Electron n: " << p.n << " Timer: " << p.exciteTimer << endl;
           }
           for (Wave& w : waves) {
               p.checkAbsorbtion(w);
           }
          
       }

       // --- Draw Waves ---
       for (Wave& w : waves) {
           if (w.energy <= 0.0f) continue;
           w.Draw();
           w.update(dt);
       }

       glfwSwapBuffers(engine.window);
       glfwPollEvents();
   }


   return 0;
}




// float useX = fmodf(p.localPos.x, 10.0f);
// useX -= 5.0f;


// float envelope = expf(-(useX * useX) / (2.0f * sigma * sigma));
// float carrier  = a * sinf(k * useX - phase);
// p.localPos.y = envelope * carrier;





