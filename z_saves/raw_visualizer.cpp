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

const float a0 = 52.9f; // Bohr radius in pm


struct Camera {
   vec3 target = vec3(0.0f);
   vec3 up = vec3(0,1,0);

   // make radius positive and on a sensible scale so zoom is meaningful
   float radius = 5000.0f;

   float azimuth = 0.0f;
   float elevation = M_PI / 2.0f;

   float orbitSpeed = 0.01f;
   double zoomFactor = 2.0;

   bool dragging = false;
   double lastX = 0.0, lastY = 0.0;

   vec3 getPosition() const {
       float x = radius * sin(elevation) * cos(azimuth);
       float y = radius * cos(elevation);
       float z = radius * sin(elevation) * sin(azimuth);
       return target + vec3(x, y, z);
   }

   void rotate(float dAzimuth, float dElevation) {
       azimuth += dAzimuth;
       elevation += dElevation;


       if(elevation < 0.01f) elevation = 0.01f;
       if(elevation > M_PI - 0.01f) elevation = M_PI - 0.01f;


       if(azimuth > 2*M_PI) azimuth -= 2*M_PI;
       if(azimuth < 0) azimuth += 2*M_PI;
   }

   void processMouseMove(double x, double y) {
       if(dragging) {
           double dx = x - lastX;
           double dy = y - lastY;
           rotate(-dx * orbitSpeed, -dy * orbitSpeed); // fixed axes
       }
       lastX = x;
       lastY = y;
   }
   void processMouseButton(int button, int action, int mods, GLFWwindow* win) {
       if(button == GLFW_MOUSE_BUTTON_LEFT) {
           if(action == GLFW_PRESS) {
               dragging = true;
               // get cursor in GLFW coords (origin top-left) then convert to OpenGL coords (origin bottom-left)
               double lx, ly;
               glfwGetCursorPos(win, &lx, &ly);
               int ww, wh;
               glfwGetWindowSize(win, &ww, &wh);
               lastX = lx;
               lastY = (double)wh - ly;
           } else if(action == GLFW_RELEASE) {
               dragging = false;
           }
       }
   }

   void processScroll(double xoffset, double yoffset) {
      // multiplicative zoom — smaller step for finer control
      const float zoomStep = 0.10f; // 10% per notch
      float factor;
      if (yoffset > 0.0) factor = powf(1.0f - zoomStep, (float) yoffset);   // zoom in
      else factor = powf(1.0f + zoomStep, (float) -yoffset);               // zoom out
      radius *= factor;


      if (radius < 1.0f) radius = 1.0f;
      if (radius > 50000.0f) radius = 50000.0f;
  }
   vec2 projectTo2D(const vec3& p, float left, float right, float bottom, float top) const {
       vec3 camPos = getPosition();
       mat4 view = lookAt(camPos, target, up);
       mat4 proj = ortho(left, right, bottom, top, -1.0f, 1.0f);


       vec4 clip = proj * view * vec4(p, 1.0f);


       // NDC
       float xN = clip.x / clip.w;
       float yN = clip.y / clip.w;


       // flip Y so result uses bottom-left origin like your drawing code
       yN = -yN;


       // Scale to your engine units (match engine.run projection below)
       float x = left + (xN + 1.0f) * 0.5f * (right - left);
       float y = bottom + (yN + 1.0f) * 0.5f * (top - bottom);


       return vec2(x, y);
   }
};
Camera camera;


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


       // Use framebuffer size for viewport to support high-DPI displays
       int fbWidth, fbHeight;
       glfwGetFramebufferSize(window, &fbWidth, &fbHeight);
       glViewport(0, 0, fbWidth, fbHeight);
   };
   void run() {
       glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


       // Use camera.radius to control orthographic extents (zoom)
       // Choose a base scale so radius==200 gives the original extents


       const float baseRadius = 200.0f;
       float zoomScale = camera.radius / baseRadius;
       if (zoomScale <= 0.001f) zoomScale = 0.001f;
       double left   = -WIDTH * zoomScale;
       double right  =  WIDTH * zoomScale;
       double bottom = -HEIGHT * zoomScale;
       double top    =  HEIGHT * zoomScale;


       // Save current matrices
       glMatrixMode(GL_PROJECTION);
       glPushMatrix();
       glLoadIdentity();
       glOrtho(left, right, bottom, top, -1.0, 1.0);
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

   vec3 sphericalToCartesian(float r, float theta, float phi){
       float x = r * sin(theta) * cos(phi);
       float y = r * cos(theta);
       float z = r * sin(theta) * sin(phi);
       return vec3(x, y, z);
   };
};
Engine engine;


struct Particle {
    vec3 pos;
    vec3 color = vec3(1.0f, 0.0f, 0.0f);

    Particle(const vec3& pos, const vec3& col) : pos(pos), color(col) {}
};
struct Particle_2d {
    vec2 pos;
    vec3 color = vec3(0.0f, 1.0f, 0.0f);

    Particle_2d(const vec2& pos, const vec3& col) : pos(pos), color(col) {}
};
struct MapVal {
    float density = 0.0f;
    vec3 color = vec3(1.0f, 1.0f, 1.0f);
    vec3 colorSum = vec3(0.0f);
};



// ================= Callbacks ================= //
void setupCameraCallbacks(GLFWwindow* window) {
   glfwSetWindowUserPointer(window, &camera);


   glfwSetMouseButtonCallback(window, [](GLFWwindow* win, int button, int action, int mods) {
       Camera* cam = (Camera*)glfwGetWindowUserPointer(win);
       cam->processMouseButton(button, action, mods, win);
   });

   // flip Y here so Camera sees bottom-left origin coordinates
   glfwSetCursorPosCallback(window, [](GLFWwindow* win, double x, double y) {
       Camera* cam = (Camera*)glfwGetWindowUserPointer(win);
       int ww, wh;
       glfwGetWindowSize(win, &ww, &wh);
       double yf = (double)wh - y;
       cam->processMouseMove(x, yf);
   });

   glfwSetScrollCallback(window, [](GLFWwindow* win, double xoffset, double yoffset) {
       Camera* cam = (Camera*)glfwGetWindowUserPointer(win);
       cam->processScroll(xoffset, yoffset);
   });
}




// ================= Grid Vars ================= //
const int   GRID_W = engine.WIDTH;
const int   GRID_H = engine.HEIGHT;
const float CELL_W = 1.0f;
const float CELL_H = 1.0f;


vector<MapVal> density(GRID_W * GRID_H);
float maxDensity = 1.0f;

vector<Particle> LoadWavefunction(const string& filename) {
    vector<Particle> pts;
    std::ifstream file("orbitals/" + filename);
    if (!file.is_open()) {
        cerr << "Failed to open JSON file: " << filename << endl;
        return pts;
    }

    json j;
    file >> j;

    const float bohr_to_pm = 52.9f;

    for (auto& point : j["points"]) {
        float x = point[0].get<float>() * bohr_to_pm;
        float y = point[1].get<float>() * bohr_to_pm;
        float z = point[2].get<float>() * bohr_to_pm;

        // Particle radius (small for electrons)
        float radius = 1.0f;

        // Color: blue for electron
        vec4 color = vec4(0.2f, 0.5f, 1.0f, 1.0f);
        //vec4 color = vec4(1.0f, 1.0f, 0.0f, 1.0f);

        pts.emplace_back(vec3(x, y, z), color);
    }

    return pts;
}

// ================= Density map functions ================= //
vec3 densityToColour(float d, float maxD, vec3 color) {

    //cout<<"color: "<<color.r<<", "<<color.g<<", "<<color.b<<endl;

    if (maxD <= 0.0f) maxD = 1.0f;
    float t = glm::clamp(d / maxD, 0.0f, 1.0f);
    float r, g, b;  
    const vec3 black(0.0f, 0.0f, 0.0f);
    const vec3 white(1.0f, 1.0f, 1.0f); 

    if (t < 0.5f) {
        // black -> cyan
        float k = t / 0.5f;          // 0..1
        return glm::mix(black, color, k);
    } else {
        // cyan -> white
        float k = (t - 0.5f) / 0.5f; // 0..1
        return glm::mix(color, white, k);
        // float k = t / 0.5f;          // 0..1
        // return glm::mix(black, cyan, k);
    }   
    return (vec3(r / 1.0, g / 1.0, b / 1.0));
}
void calculateDensityMap(const vector<Particle_2d>& particles) {
   // Reset
   for (MapVal& val : density) {
       val.density = 0.0f;
   }
   maxDensity = 0.0f;


   // Match ortho extents used by engine.run
   const float baseRadius = 200.0f;
   float zoomScale = camera.radius / baseRadius;
   if (zoomScale <= 0.001f) zoomScale = 0.001f;
   float left   = -engine.WIDTH * zoomScale;
   float right  =  engine.WIDTH * zoomScale;
   float bottom = -engine.HEIGHT * zoomScale;
   float top    =  engine.HEIGHT * zoomScale;


   // Cell size in world coordinates
   float cellW = (right - left) / float(GRID_W);
   float cellH = (top - bottom) / float(GRID_H);


   // Influence radius in world units (tweak if you want wider/narrower blur)
   const float INFLUENCE_RADIUS = 200.0f;
   const float INFLUENCE_RADIUS_SQ = INFLUENCE_RADIUS * INFLUENCE_RADIUS;


   for (const Particle_2d& par : particles) {
        vec2 p = par.pos;
        // particle outside current extents -> skip
        if (p.x < left || p.x > right || p.y < bottom || p.y > top) continue;


        int center_i = int(floor((p.x - left) / cellW));
        int center_j = int(floor((p.y - bottom) / cellH));
        if (center_i < 0 || center_i >= GRID_W || center_j < 0 || center_j >= GRID_H) continue;
        int radiusCellsX = (int)ceil(INFLUENCE_RADIUS / cellW);
        int radiusCellsY = (int)ceil(INFLUENCE_RADIUS / cellH);
        int min_i = std::max(0, center_i - radiusCellsX);
        int max_i = std::min(GRID_W - 1, center_i + radiusCellsX);
        int min_j = std::max(0, center_j - radiusCellsY);
        int max_j = std::min(GRID_H - 1, center_j + radiusCellsY);
        for (int j = min_j; j <= max_j; ++j) {
            for (int i = min_i; i <= max_i; ++i) {
                // cell center in world coords
                float cell_x = left + (i + 0.5f) * cellW;
                float cell_y = bottom + (j + 0.5f) * cellH;
                float dx = p.x - cell_x;
                float dy = p.y - cell_y;
                float distSq = dx*dx + dy*dy;
                if (distSq < INFLUENCE_RADIUS_SQ) {
                    float dist = sqrtf(distSq);
                    float influence = (INFLUENCE_RADIUS - dist) / INFLUENCE_RADIUS;
                    float kernel_weight = influence * influence;
                    int index = j * GRID_W + i;
                    density[index].density += kernel_weight;
                    density[index].colorSum += par.color * 0.004f;
                    // CRITICAL FIX: Accumulate weighted color
                    //density[index].colorSum += par.color * kernel_weight;
                    maxDensity = std::max(maxDensity, density[index].density);
                }
            }
        }
    }
}
void drawDensityMap() {
   // Match ortho extents used by engine.run
   const float baseRadius = 200.0f;
   float zoomScale = camera.radius / baseRadius;
   if (zoomScale <= 0.001f) zoomScale = 0.001f;
   float left   = -engine.WIDTH * zoomScale;
   float right  =  engine.WIDTH * zoomScale;
   float bottom = -engine.HEIGHT * zoomScale;
   float top    =  engine.HEIGHT * zoomScale;


   float cellW = (right - left) / float(GRID_W);
   float cellH = (top - bottom) / float(GRID_H);


   glBegin(GL_QUADS);
   for (int j = 0; j < GRID_H; ++j) {
       for (int i = 0; i < GRID_W; ++i) {
           int index = j * GRID_W + i;
           float d = density[index].density;

            
           vec3 color = densityToColour(d, maxDensity, density[index].colorSum);
           glColor3f(color.r, color.g, color.b);


           float x0 = left + i * cellW;
           float y0 = bottom + j * cellH;
           float x1 = x0 + cellW;
           float y1 = y0 + cellH;


           glVertex2f(x0, y0);
           glVertex2f(x1, y0);
           glVertex2f(x1, y1);
           glVertex2f(x0, y1);
       }
   }
   glEnd();
}


void project_2d(const vector<Particle>& particles_3d, vector<Particle_2d>& particles_2d) {
   particles_2d.clear();


   // match the ortho extents used in Engine::run so projection/viewport align
   const float baseRadius = 200.0f;
   float zoomScale = camera.radius / baseRadius;
   if (zoomScale <= 0.001f) zoomScale = 0.001f;


   float left   = -engine.WIDTH * zoomScale;
   float right  =  engine.WIDTH * zoomScale;
   float bottom = -engine.HEIGHT * zoomScale;
   float top    =  engine.HEIGHT * zoomScale;


   for (const Particle& p : particles_3d) {
        particles_2d.push_back( 
            Particle_2d(camera.projectTo2D(p.pos, left, right, bottom, top) , p.color)
        );
   }
}
void generateParticles(vector<vec3>& particles, int numParticles, int numClusters, float clusterRadius)
{
   particles.clear();


   // 1. Simulation bounds converted to 3D
   float minX = -engine.WIDTH;
   float maxX = engine.WIDTH;


   float minY = -engine.HEIGHT;
   float maxY = engine.HEIGHT;


   float minZ = -500;     // add this to your engine
   float maxZ = 500;


   // 2. Randomly pick cluster centers inside the 3D bounding box
   std::vector<vec3> clusterCenters;
   clusterCenters.reserve(numClusters);


   for (int i = 0; i < numClusters; ++i) {
       float x = minX + (std::rand() / (float)RAND_MAX) * (maxX - minX);
       float y = minY + (std::rand() / (float)RAND_MAX) * (maxY - minY);
       float z = minZ + (std::rand() / (float)RAND_MAX) * (maxZ - minZ);
       clusterCenters.emplace_back(x, y, z);
   }


   // 3. Global (uniform) vs clustered particles
   int globalParticles = numParticles / 3;
   int clusteredParticles = numParticles - globalParticles;


   // --- Global random 3D particles ---
   for (int i = 0; i < globalParticles; ++i) {
       float x = minX + (std::rand() / (float)RAND_MAX) * (maxX - minX);
       float y = minY + (std::rand() / (float)RAND_MAX) * (maxY - minY);
       float z = minZ + (std::rand() / (float)RAND_MAX) * (maxZ - minZ);
       particles.emplace_back(x, y, z);
   }


   // --- Clustered 3D particles (uniform spherical distribution) ---
   for (int i = 0; i < clusteredParticles; ++i) {
       const vec3& center = clusterCenters[std::rand() % numClusters];


       // Random direction on a sphere
       float u = (std::rand() / (float)RAND_MAX) * 2.0f - 1.0f; // cos(theta)
       float phi = (std::rand() / (float)RAND_MAX) * 2.0f * M_PI;


       float sqrt1MinusU2 = sqrt(1 - u * u);
       float dx = sqrt1MinusU2 * cos(phi);
       float dy = sqrt1MinusU2 * sin(phi);
       float dz = u;


       // Random radius inside sphere (cube-root preserves uniform density)
       float radius = clusterRadius * cbrt(std::rand() / (float)RAND_MAX);


       particles.emplace_back(
           center.x + dx * radius,
           center.y + dy * radius,
           center.z + dz * radius
       );
   }
}
void initRandom() {
   srand(static_cast<unsigned int>(time(0)));
}

// Generate a random unit vector
vec3 randomUnitVector() {
   float theta = static_cast<float>(rand()) / RAND_MAX * 2.0f * glm::pi<float>(); // 0 to 2pi
   float phi   = acos(2.0f * static_cast<float>(rand()) / RAND_MAX - 1.0f);       // 0 to pi
   float x = sin(phi) * cos(theta);
   float y = sin(phi) * sin(theta);
   float z = cos(phi);
   return vec3(x, y, z);
}
void generateDirs(vec3& dir1, vec3& dir2, float angleDeg) {
   dir1 = randomUnitVector();


   // Rotate dir1 by angleDeg around a random perpendicular axis
   vec3 axis = normalize(cross(dir1, randomUnitVector()));
   float angleRad = glm::radians(angleDeg);
   dir2 = dir1 * cos(angleRad) + cross(axis, dir1) * sin(angleRad);
   // dir2 is automatically unit length since both dir1 and axis are unit
}
void generateTetrahedralDirs(vec3 dirs[4]) {
   float a = sqrt(3.0f) / 3.0f; // magnitude factor for tetrahedron corners
   dirs[0] = normalize(vec3( a,  a,  a));
   dirs[1] = normalize(vec3(-a, -a,  a));
   dirs[2] = normalize(vec3(-a,  a, -a));
   dirs[3] = normalize(vec3( a, -a, -a));
}


// ================= Main ================= //
int main () {
    setupCameraCallbacks(engine.window);
    
    vector<Particle> particles_3d = LoadWavefunction("orbital_n3_l2_m0.json");

    // ------- 1. Declare 2D projection -------------------
    vector<Particle_2d> particles;


    while (!glfwWindowShouldClose(engine.window)) {
        engine.run();

        // -------- project 3D to 2D and recalculate density --------
        project_2d(particles_3d, particles);

        for (MapVal & val : density){
            val.density = 0.0f;
            val.colorSum = vec3(0.0f);
        }

        //calculateDensityMap( particles );
        //calculateDensityMap(particles3p_x_2d, density3p_x);
       


        // // ------- 1. Draw the density map -------------------
        //drawDensityMap();

        //break;
        // ------- 2. Draw the particles -------------------
        // Draw 2D projected particles in this orthographic space
        glColor4f(1.0f, 1.0f, 1.0f, 0.7f);
       
        for (const Particle_2d& p : particles){
            engine.drawFilledCircle(p.pos.x, p.pos.y, 1.0f);
        }


        glfwSwapBuffers(engine.window);
        glfwPollEvents();
    }


    return 0;
}



