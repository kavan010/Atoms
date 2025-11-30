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


struct Atom {
    
    
    vec3 pos;
    int atomicNum = 0;
    
    
    Atom(vec3 pos, int atomicNum) : pos(pos), atomicNum(atomicNum) {}
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
    float radialProbability3s(float r) {
        float x = r / a0;
        float R = (1.0f - (2.0f/3.0f)*x + (2.0f/27.0f)*x*x) * exp(-x / 3.0f);
        return R * R * r * r;
    }
    float radialProbability3p(float r) {
        float x = r / a0;
        float R = x * (1.0f - x / 6.0f) * exp(-x / 3.0f);
        return R * R * r * r; // multiply by r^2 for probability density
    }
    float radialProbability3dxy(float r) {
        float x = r / a0;
        float R = x*x * (6.0f - x) * exp(-x / 3.0f); // 3d radial part
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
    float sampleR3s() {
        float r_max = 20.0f * a0;  // 3s extends farther than 1s/2s
        float P_max = radialProbability3s(3.0f * a0); // near the 3s radial peak
    
    
        while (true) {
            float r = static_cast<float>(rand()) / RAND_MAX * r_max;
            float y = static_cast<float>(rand()) / RAND_MAX * P_max;
            if (y <= radialProbability3s(r)) return r;
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
    float sampleR3d() {
        // set a generous r_max
        float r_max = 30.0f * a0;
        // true peak for 3d radial probability occurs at r = 9 a0
        float P_max = radialProbability3dxy(9.0f * a0);

        while (true) {
            float u1 = static_cast<float>(rand()) / RAND_MAX; // [0,1)
            float r = u1 * r_max;
            float u2 = static_cast<float>(rand()) / RAND_MAX;
            float y = u2 * P_max;
            if (y <= radialProbability3dxy(r)) return r;
        }
    }


   void sample1s(vector<Particle> &particles) {   // change return type to void
       float r = sampleR1s();
       float theta = acos(1.0f - 2.0f * static_cast<float>(rand()) / RAND_MAX); // [0, pi]
       float phi   = 2.0f * M_PI * static_cast<float>(rand()) / RAND_MAX;       // [0, 2pi]
       vec3 electronPos = engine.sphericalToCartesian(r, theta, phi) + pos;
       vec3 color(0.0f, 1.0f, 1.0f); // cyan-ish color for 1s
      
       // Construct Particle in-place
       if (true) {
           particles.emplace_back(  Particle(electronPos, color ) ); // cyan-ish color
       }
   }
   void sample2s(vector<Particle> &particles) {
       float r = sampleR2s();
       float theta = acos(1.0f - 2.0f * static_cast<float>(rand()) / RAND_MAX); // [0, π]
       float phi   = 2.0f * M_PI * static_cast<float>(rand()) / RAND_MAX;       // [0, 2π]
      
       vec3 electronPos = engine.sphericalToCartesian(r, theta, phi) + pos;
       vec3 color(1.0f, 0.0f, 1.0f); // magenta-ish color

       // Use a distinct color for visualization, e.g. yellow for 2s
       if (true) {
       particles.emplace_back( Particle(electronPos, color )); // cyan-ish color
       }
   }
   void sample3s(vector<Particle> &particles) {
        float r = sampleR3s();

        // Uniform spherical angle distribution
        float theta = acos(1.0f - 2.0f * static_cast<float>(rand()) / RAND_MAX);
        float phi   = 2.0f * M_PI * static_cast<float>(rand()) / RAND_MAX;

        vec3 pos = engine.sphericalToCartesian(r, theta, phi) + pos;
        vec3 color(1.0f, 1.0f, 0.0f); // yellow-ish color for 3s

        particles.emplace_back( Particle(pos, color) );
   }
   void sample2p_x(vector<Particle> &particles) {
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


       vec3 electronPos = engine.sphericalToCartesian(r, theta, phi) + pos;
       vec3 color(1.0f, 0.0f, 0.0f); // red color for 2p_x


       // Keep one lobe for visualization
       if (true) {
           particles.emplace_back(  Particle(electronPos, color) ); // red for 2p_x
       }
   }
   void sample2p_y(vector<Particle> &particles) {
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


       vec3 electronPos = engine.sphericalToCartesian(r, theta, phi) + pos;
       vec3 color(1.0f, 0.0f, 0.0f); // red color for 2p_y


       // Keep one lobe for visualization
       if (true) {
           particles.emplace_back( Particle(electronPos, color) ); // red for 2p_y (keeps existing color)
       }
   }
   void sample2p_z(vector<Particle> &particles) {
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


       vec3 electronPos = engine.sphericalToCartesian(r, theta, phi) + pos;
       vec3 color(1.0f, 0.0f, 0.0f); // red color for 2p_z

       // Keep one lobe for visualization
       particles.emplace_back(Particle(electronPos, color)); // red for 2p_z
   }

   void sample3p_x(vector<Particle> &particles) {
       float r = sampleR3p(); // 3p radial distribution
       float theta, phi;


       // sample theta and phi
       while (true) {
           theta = acos(1.0f - 2.0f * static_cast<float>(rand()) / RAND_MAX); // [0, pi]
           phi   = 2.0f * M_PI * static_cast<float>(rand()) / RAND_MAX;       // [0, 2pi]


           float prob = pow(sin(theta) * cos(phi), 2); // sin^2(theta) * cos^2(phi)
           if (static_cast<float>(rand()) / RAND_MAX <= prob) break;
       }


       vec3 electronPos = engine.sphericalToCartesian(r, theta, phi) + pos;
       vec3 color(1.0f, 0.0f, 1.0f); // blue-ish

       particles.emplace_back(Particle(electronPos, color)); // assign a color for 3p_x
   }
   void sample3p_y(vector<Particle> &particles) {
       float r = sampleR3p(); // 3p radial distribution
       float theta, phi;


       // sample theta and phi
       while (true) {
           theta = acos(1.0f - 2.0f * static_cast<float>(rand()) / RAND_MAX); // [0, pi]
           phi   = 2.0f * M_PI * static_cast<float>(rand()) / RAND_MAX;       // [0, 2pi]


           float prob = pow(sin(theta) * sin(phi), 2); // sin^2(theta) * sin^2(phi)
           if (static_cast<float>(rand()) / RAND_MAX <= prob) break;
       }


       vec3 electronPos = engine.sphericalToCartesian(r, theta, phi) + pos;
       vec3 color(0.7f, 0.0f, 1.0f); // blue-ish

       particles.emplace_back(Particle(electronPos, color)); // assign a color for 3p_y
   }
   void sample3p_z(vector<Particle> &particles) {
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


       vec3 electronPos = engine.sphericalToCartesian(r, theta, phi) + pos;
       vec3 color(1.0f, 1.0f, 0.0f);


       // visualize (different color for 3p_z)
       particles.emplace_back( Particle(electronPos, color) ); // cyan-ish for 3p_z
   }
  
    void sample3dxy( vector<Particle> &particles) {
        float r = sampleR3d();

        float theta, phi;
        while (true) {
            float u = static_cast<float>(rand()) / RAND_MAX;
            float v = static_cast<float>(rand()) / RAND_MAX;
            theta = acosf(1.0f - 2.0f * u);
            phi = 2.0f * M_PI * v; 
            
            float s = sinf(theta);
            float angProb = (s * s) * (s * s) * (sinf(2.0f * phi) * sinf(2.0f * phi));

            float rcheck = static_cast<float>(rand()) / RAND_MAX;
            if (rcheck <= angProb) break;
        }

        vec3 electronPos = engine.sphericalToCartesian(r, theta, phi);

        vec3 color(0.0f, 1.0f, 1.0f);
        particles.emplace_back( Particle(electronPos, color) ); 
    }
    void sample3dxz( vector<Particle> &particles) {
        float r = sampleR3d();

        float theta, phi;
        while (true) {
            float u = float(rand()) / RAND_MAX;
            float v = float(rand()) / RAND_MAX;

            theta = acosf(1.0f - 2.0f * u);
            phi = 2.0f * M_PI * v;

            float s = sinf(theta);
            float c = cosf(theta);

            // Angular probability for d_xz:
            // |ψ|² ∝ sin²θ * cos²θ * cos²φ
            float angProb =
                (s * s) * (c * c) *
                (cosf(phi) * cosf(phi));

            if ((float(rand()) / RAND_MAX) <= angProb)
                break;
        }

        vec3 electronPos = engine.sphericalToCartesian(r, theta, phi);
        vec3 color(0.0f, 1.0f, 1.0f);
        particles.emplace_back( Particle(electronPos, color) ); 
    }
    void sample3dyz( vector<Particle> &particles) {
        float r = sampleR3d();

        float theta, phi;
        while (true) {
            float u = float(rand()) / RAND_MAX;
            float v = float(rand()) / RAND_MAX;

            theta = acosf(1.0f - 2.0f * u);
            phi = 2.0f * M_PI * v;

            float s = sinf(theta);
            float c = cosf(theta);

            // Angular probability for d_yz:
            // |ψ|² ∝ sin²θ * cos²θ * sin²φ
            float angProb =
                (s * s) * (c * c) *
                (sinf(phi) * sinf(phi));

            if ((float(rand()) / RAND_MAX) <= angProb)
                break;
        }

        vec3 electronPos = engine.sphericalToCartesian(r, theta, phi);
        vec3 color(0.0f, 1.0f, 1.0f);
        particles.emplace_back( Particle(electronPos, color) ); 
    }

};
vector<Atom> atoms = {

   //Atom(vec3(0.0f, 0.0f, 0.0f), 6),
   // Atom(vec3(0.0f, 2500.0f, 0.0f))
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


// ======= molecules =======
void addCO2(vec3 centre) {
   atoms.emplace_back(centre, 6); // Carbon


   float bondLength = 465.0f;
   vec3 dir1, dir2;
   generateDirs(dir1, dir2, 180.0f); // CO2 is linear, O-C-O ~180°


   atoms.emplace_back(centre + dir1 * bondLength, 8); // Oxygen
   atoms.emplace_back(centre + dir2 * bondLength, 8); // Oxygen
}
void addCH4(vec3 centre) {
   atoms.emplace_back(centre, 6); // Carbon


   float bondLength = 410.0f; // approximate C-H bond length in your units
   vec3 dirs[4];


   // Generate four directions for tetrahedral geometry (~109.5° between bonds)
   generateTetrahedralDirs(dirs);


   // Place Hydrogens
   for (int i = 0; i < 4; ++i) {
       atoms.emplace_back(centre + dirs[i] * bondLength, 1); // Hydrogen
   }
}
void addWater(vec3 centre) {
   atoms.emplace_back(centre, 8); // Oxygen


   float bondLength = 95.9f * 4.5f; // H-O bond ~0.96Å scaled up


   vec3 dir1, dir2;
   generateDirs(dir1, dir2, 104.5f); // H-O-H bond ~104.5°


   atoms.emplace_back( centre + dir1 * bondLength, 1); // Hydrogen
   atoms.emplace_back(centre + dir2 * bondLength, 1); // Hydrogen
}


// ================= Main ================= //
int main () {
    setupCameraCallbacks(engine.window);
    
    const int numMolecules = 50;   // number of water molecules
    const float spacing = 2000.0f; // approximate distance between molecules
    vector<vec3> molecules;


    //addCH4(vec3(0.0f));
    // addCH4(vec3(0.0f, 0.0f, 0.0f));
    atoms = {
        Atom(vec3(0.0f), 102)
    };
    //addWater(vec3(0.0f, 0.0f, 2000.0f));
    //addCO2(vec3(0.0f, 0.0f, -2000.0f));
    

    // ----- Generate particles -----
    vector<Particle> particles_3d;
    vector<Particle> particles3p_x_3d;
    //generateParticles(particles_3d, 10000, 5, 300.0f);
    for (Atom & atom : atoms) {
        for (int i = 0; i < 5000; i++) {
            if (atom.atomicNum == 1) 
                atom.sample1s(particles_3d);
            else if (atom.atomicNum == 6) {
                atom.sample1s(particles_3d);
                atom.sample1s(particles_3d);
                atom.sample2s(particles_3d);
                atom.sample2s(particles_3d);
                atom.sample2p_x(particles_3d);
                atom.sample2p_y(particles_3d);
            } else if (atom.atomicNum == 8) {
                atom.sample1s(particles_3d);
                atom.sample1s(particles_3d);
                atom.sample2s(particles_3d);
                atom.sample2s(particles_3d);
                atom.sample2p_x(particles_3d);
                atom.sample2p_x(particles_3d);
                atom.sample2p_y(particles_3d);
                atom.sample2p_z(particles_3d);
            } else if(atom.atomicNum == 102)
                //atom.sample3dxy(particles_3d);
                //atom.sample3dxz(particles_3d);
                atom.sample3dyz(particles_3d);
        }
    }


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

        calculateDensityMap( particles );
        //calculateDensityMap(particles3p_x_2d, density3p_x);
       


        // // ------- 1. Draw the density map -------------------
        drawDensityMap();

        //break;
        // ------- 2. Draw the particles -------------------
        // Draw 2D projected particles in this orthographic space
        // glColor4f(1.0f, 1.0f, 1.0f, 0.7f);
       
        // for (const Particle_2d& p : particles){
        //     engine.drawFilledCircle(p.pos.x, p.pos.y, 1.0f);
        // }


        glfwSwapBuffers(engine.window);
        glfwPollEvents();
    }


    return 0;
}



