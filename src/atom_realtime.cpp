#include <GL/glew.h>
#include <GLFW/glfw3.h>
#ifndef __APPLE__
#include <GL/glu.h>
#endif
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
#include <complex>
#include <random>
#include <string>
#include <array>
#include <algorithm>
#include <cctype>
#include <unordered_map>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
using namespace glm;
using namespace std;

// ================= Constants ================= //
const float a0 = 1;
float electron_r = 1.5f; // radius for spheres
const double hbar = 1;
const double m_e = 1;
const double zmSpeed = 10.0;

// --- Global quantum numbers ---
int n = 2, l = 1, m = 0, N = 100000;

struct QuantumCombo {
    int n;
    int l;
    int m;
};

vector<QuantumCombo> quantumCombos;
int currentComboIndex = -1;
constexpr int selectorMaxN = 7;

struct ElementQuantum {
    int atomicNumber;
    string name;
    int n;
    int l;
    int m;
};

vector<ElementQuantum> elementMappings;
int currentElementIndex = -1;

int magneticQuantumForLastElectron(int lValue, int electronInSubshell) {
    const int orbitals = 2 * lValue + 1;
    if (electronInSubshell <= orbitals) {
        return -lValue + (electronInSubshell - 1);
    }

    return -lValue + (electronInSubshell - orbitals - 1);
}

QuantumCombo quantumFromAtomicNumber(int atomicNumber) {
    static const vector<pair<int, int>> orbitalOrder = {
        {1, 0}, {2, 0}, {2, 1}, {3, 0}, {3, 1}, {4, 0}, {3, 2}, {4, 1}, {5, 0}, {4, 2},
        {5, 1}, {6, 0}, {4, 3}, {5, 2}, {6, 1}, {7, 0}, {5, 3}, {6, 2}, {7, 1}, {8, 0}
    };

    int remainingElectrons = atomicNumber;

    for (const auto& orbital : orbitalOrder) {
        const int nValue = orbital.first;
        const int lValue = orbital.second;
        const int capacity = 2 * (2 * lValue + 1);
        const int usedHere = std::min(remainingElectrons, capacity);

        if (remainingElectrons <= capacity) {
            return {
                nValue,
                lValue,
                magneticQuantumForLastElectron(lValue, usedHere)
            };
        }

        remainingElectrons -= usedHere;
    }

    return {1, 0, 0};
}

void buildElementMappings() {
    static const array<const char*, 118> elementNames = {{
        "Hydrogen", "Helium", "Lithium", "Beryllium", "Boron", "Carbon", "Nitrogen", "Oxygen", "Fluorine", "Neon",
        "Sodium", "Magnesium", "Aluminum", "Silicon", "Phosphorus", "Sulfur", "Chlorine", "Argon", "Potassium", "Calcium",
        "Scandium", "Titanium", "Vanadium", "Chromium", "Manganese", "Iron", "Cobalt", "Nickel", "Copper", "Zinc",
        "Gallium", "Germanium", "Arsenic", "Selenium", "Bromine", "Krypton", "Rubidium", "Strontium", "Yttrium", "Zirconium",
        "Niobium", "Molybdenum", "Technetium", "Ruthenium", "Rhodium", "Palladium", "Silver", "Cadmium", "Indium", "Tin",
        "Antimony", "Tellurium", "Iodine", "Xenon", "Cesium", "Barium", "Lanthanum", "Cerium", "Praseodymium", "Neodymium",
        "Promethium", "Samarium", "Europium", "Gadolinium", "Terbium", "Dysprosium", "Holmium", "Erbium", "Thulium", "Ytterbium",
        "Lutetium", "Hafnium", "Tantalum", "Tungsten", "Rhenium", "Osmium", "Iridium", "Platinum", "Gold", "Mercury",
        "Thallium", "Lead", "Bismuth", "Polonium", "Astatine", "Radon", "Francium", "Radium", "Actinium", "Thorium",
        "Protactinium", "Uranium", "Neptunium", "Plutonium", "Americium", "Curium", "Berkelium", "Californium", "Einsteinium", "Fermium",
        "Mendelevium", "Nobelium", "Lawrencium", "Rutherfordium", "Dubnium", "Seaborgium", "Bohrium", "Hassium", "Meitnerium", "Darmstadtium",
        "Roentgenium", "Copernicium", "Nihonium", "Flerovium", "Moscovium", "Livermorium", "Tennessine", "Oganesson"
    }};

    elementMappings.clear();
    elementMappings.reserve(elementNames.size());

    for (size_t i = 0; i < elementNames.size(); ++i) {
        const int atomicNumber = static_cast<int>(i) + 1;
        const QuantumCombo combo = quantumFromAtomicNumber(atomicNumber);
        elementMappings.push_back({
            atomicNumber,
            elementNames[i],
            combo.n,
            combo.l,
            combo.m
        });
    }
}

int findElementIndexByQuantum(int nValue, int lValue, int mValue) {
    for (size_t i = 0; i < elementMappings.size(); ++i) {
        const ElementQuantum& element = elementMappings[i];
        if (element.n == nValue && element.l == lValue && element.m == mValue) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

void syncElementSelectionFromQuantum() {
    currentElementIndex = findElementIndexByQuantum(n, l, m);
}

void setElementByIndex(int elementIndex) {
    if (elementMappings.empty()) return;
    const int total = static_cast<int>(elementMappings.size());

    currentElementIndex = elementIndex % total;
    if (currentElementIndex < 0) currentElementIndex += total;

    const ElementQuantum& element = elementMappings[currentElementIndex];
    n = element.n;
    l = element.l;
    m = element.m;
}

void stepElementSelection(int direction) {
    if (elementMappings.empty()) return;

    if (currentElementIndex < 0) {
        syncElementSelectionFromQuantum();
        if (currentElementIndex < 0) currentElementIndex = 0;
    }

    setElementByIndex(currentElementIndex + direction);
}

string currentElementLabel() {
    if (currentElementIndex >= 0 && currentElementIndex < static_cast<int>(elementMappings.size())) {
        const ElementQuantum& element = elementMappings[currentElementIndex];
        return to_string(element.atomicNumber) + " " + element.name;
    }

    return "Custom Orbital";
}

// ================= Physics Sampling ================= //
struct Particle {
    vec3 pos;
    vec3 vel = vec3(0.0f);
    vec4 color;
    Particle(vec3 p, vec4 c = vec4(0.0f, 0.5f, 1.0f, 1.0f)) : pos(p), color(c){}
};
vector<Particle> particles;

// --- random devices ---
random_device rd; mt19937 gen(rd()); uniform_real_distribution<float> dis(0.0f, 1.0f);

// --- sample R ---                 <- uses CDF sampling
double sampleR(int n, int l, mt19937& gen) {
    const int N = 4096;
    //const double a0 = 1.0;
    const double rMax = 10.0 * n * n * a0;

    static vector<double> cdf;
    static bool built = false;
    static int cachedN = -1;
    static int cachedL = -1;

    if (!built || cachedN != n || cachedL != l) {
        cdf.resize(N);
        double dr = rMax / (N - 1);
        double sum = 0.0;

        for (int i = 0; i < N; ++i) {
            double r = i * dr;
            double rho = 2.0 * r / (n * a0);

            // Associated Laguerre L_{n-l-1}^{2l+1}(rho)
            int k = n - l - 1;
            int alpha = 2 * l + 1;

            double L = 1.0, Lm1 = 1.0 + alpha - rho;
            if (k == 1) L = Lm1;
            else if (k > 1) {
                double Lm2 = 1.0;
                for (int j = 2; j <= k; ++j) {
                    L = ((2*j - 1 + alpha - rho) * Lm1 -
                         (j - 1 + alpha) * Lm2) / j;
                    Lm2 = Lm1;
                    Lm1 = L;
                }
            }

            double norm = pow(2.0 / (n * a0), 3) * tgamma(n - l) / (2.0 * n * tgamma(n + l + 1));
            double R = sqrt(norm) * exp(-rho / 2.0) * pow(rho, l) * L;

            double pdf = r * r * R * R;
            sum += pdf;
            cdf[i] = sum;
        }

        for (double& v : cdf) v /= sum;
        built = true;
        cachedN = n;
        cachedL = l;
    }

    uniform_real_distribution<double> dis(0.0, 1.0);
    double u = dis(gen);

    int idx = lower_bound(cdf.begin(), cdf.end(), u) - cdf.begin();
    return idx * (rMax / (N - 1));
}
// --- sample Theta ---             <- uses CDF sampling
double sampleTheta(int l, int m, mt19937& gen) {
    const int N = 2048;
    static vector<double> cdf;
    static bool built = false;
    static int cachedL = -1;
    static int cachedAbsM = -1;
    const int absM = abs(m);

    if (!built || cachedL != l || cachedAbsM != absM) {
        cdf.resize(N);
        double dtheta = M_PI / (N - 1);
        double sum = 0.0;

        for (int i = 0; i < N; ++i) {
            double theta = i * dtheta;
            double x = cos(theta);

            // Associated Legendre P_l^m(x)
            double Pmm = 1.0;
            if (absM > 0) {
                double somx2 = sqrt((1.0 - x) * (1.0 + x));
                double fact = 1.0;
                for (int j = 1; j <= absM; ++j) {
                    Pmm *= -fact * somx2;
                    fact += 2.0;
                }
            }

            double Plm;
            if (l == absM) {
                Plm = Pmm;
            } else {
                double Pm1m = x * (2 * absM + 1) * Pmm;
                if (l == absM + 1) {
                    Plm = Pm1m;
                } else {
                    double Pll;
                    for (int ll = absM + 2; ll <= l; ++ll) {
                        Pll = ((2 * ll - 1) * x * Pm1m -
                               (ll + absM - 1) * Pmm) / (ll - absM);
                        Pmm = Pm1m;
                        Pm1m = Pll;
                    }
                    Plm = Pm1m;
                }
            }

            double pdf = sin(theta) * Plm * Plm;
            sum += pdf;
            cdf[i] = sum;
        }

        for (double& v : cdf) v /= sum;
        built = true;
        cachedL = l;
        cachedAbsM = absM;
    }

    uniform_real_distribution<double> dis(0.0, 1.0);
    double u = dis(gen);

    int idx = lower_bound(cdf.begin(), cdf.end(), u) - cdf.begin();
    return idx * (M_PI / (N - 1));
}
// --- sample Phi (uniform) ---     <- uses CDF sampling
float samplePhi(float n, float l, float m) {
    return 2.0f * M_PI * dis(gen);
}
// --- calculate prob current ---
vec3 calculateProbabilityFlow(Particle& p, int n, int l, int m) {
    double r = length(p.pos);   if (r < 1e-6) return vec3(0.0f);
    double theta = acos(p.pos.y / r); 
    double phi = atan2(p.pos.z, p.pos.x); 


    //Compute magnitude
    double sinTheta = sin(theta);  if (abs(sinTheta) < 1e-4) sinTheta = 1e-4;
    double v_mag = hbar * m / (m_e * r * sinTheta);

    //Convert to Cartesian
    double vx = -v_mag * sin(phi);
    double vy = 0.0; 
    double vz =  v_mag * cos(phi);

    return vec3((float)vx, (float)vy, (float)vz);
}

// --- color map ---
vec4 heatmap_fire(float value) {
    // Ensure value is clamped between 0 and 1
    value = std::max(0.0f, std::min(1.0f, value));

    // Define color stops for the "Heat/Fire" pattern
    // Order: Black -> Dark Purple -> Red -> Orange -> Yellow -> White
    const int num_stops = 6;
    vec4 colors[num_stops] = {
        {0.0f, 0.0f, 0.0f, 1.0f}, // 0.0: Black
        {0.5f, 0.0f, 0.99f, 1.0f}, // 0.2: Dark Purple
        {0.8f, 0.0f, 0.0f, 1.0f}, // 0.4: Deep Red
        {1.0f, 0.5f, 0.0f, 1.0f}, // 0.6: Orange
        {1.0f, 1.0f, 0.0f, 1.0f}, // 0.8: Yellow
        {1.0f, 1.0f, 1.0f, 1.0f}  // 1.0: White
    };

    // Find which segment the value falls into
    float scaled_v = value * (num_stops - 1);
    int i = static_cast<int>(scaled_v);
    int next_i = std::min(i + 1, num_stops - 1);
    
    // Calculate how far we are between stop 'i' and 'next_i'
    float local_t = scaled_v - i;

    // Linearly interpolate between the two colors
    vec4 result;
    result.r = colors[i].r + local_t * (colors[next_i].r - colors[i].r);
    result.g = colors[i].g + local_t * (colors[next_i].g - colors[i].g);
    result.b = colors[i].b + local_t * (colors[next_i].b - colors[i].b);
    result.a = 1.0f; // Solid opacity

    return result;
}
vec4 inferno(double r, double theta, double phi, int n, int l, int m) {
    // --- radial part |R(r)|^2 ---
    double rho = 2.0 * r / (n * a0);

    int k = n - l - 1;
    int alpha = 2 * l + 1;

    double L = 1.0;
    if (k == 1) {
        L = 1.0 + alpha - rho;
    } else if (k > 1) {
        double Lm2 = 1.0;
        double Lm1 = 1.0 + alpha - rho;
        for (int j = 2; j <= k; ++j) {
            L = ((2*j - 1 + alpha - rho) * Lm1 -
                 (j - 1 + alpha) * Lm2) / j;
            Lm2 = Lm1;
            Lm1 = L;
        }
    }

    double norm = pow(2.0 / (n * a0), 3)
                * tgamma(n - l)
                / (2.0 * n * tgamma(n + l + 1));

    double R = sqrt(norm) * exp(-rho / 2.0) * pow(rho, l) * L;
    double radial = R * R;

    // --- angular part |P_l^m(cosθ)|^2 ---
    double x = cos(theta);
    const int absM = abs(m);

    double Pmm = 1.0;
    if (absM > 0) {
        double somx2 = sqrt((1.0 - x) * (1.0 + x));
        double fact = 1.0;
        for (int j = 1; j <= absM; ++j) {
            Pmm *= -fact * somx2;
            fact += 2.0;
        }
    }

    double Plm;
    if (l == absM) {
        Plm = Pmm;
    } else {
        double Pm1m = x * (2*absM + 1) * Pmm;
        if (l == absM + 1) {
            Plm = Pm1m;
        } else {
            for (int ll = absM + 2; ll <= l; ++ll) {
                double Pll = ((2*ll - 1) * x * Pm1m -
                              (ll + absM - 1) * Pmm) / (ll - absM);
                Pmm = Pm1m;
                Pm1m = Pll;
            }
            Plm = Pm1m;
        }
    }

    double angular = Plm * Plm;

    double intensity = radial * angular;

    //cout << "intensity: " << intensity << endl;

    return heatmap_fire(intensity * 1.5 * pow(5, n)); // Scale for better color mapping
}


// ================= camera ================= //
struct Camera {
    vec3 target = vec3(0.0f, 0.0f, 0.0f);
    float radius = 50.0f;
    float azimuth = 0.0f;
    float elevation = M_PI / 2.0f;
    float orbitSpeed = 0.01f;
    float panSpeed = 0.01f;
    double zoomSpeed = zmSpeed;
    bool dragging = false;
    bool panning = false;
    double lastX = 0.0, lastY = 0.0;

    vec3 position() const {
        float clampedElevation = glm::clamp(elevation, 0.01f, float(M_PI) - 0.01f);
        return vec3(
            radius * sin(clampedElevation) * cos(azimuth),
            radius * cos(clampedElevation),
            radius * sin(clampedElevation) * sin(azimuth)
        );
    }
    void update() {
        target = vec3(0.0f, 0.0f, 0.0f);
    }

    void processMouseMove(double x, double y) {
        float dx = float(x - lastX);
        float dy = float(y - lastY);
        if (dragging) {
            azimuth += dx * orbitSpeed;
            elevation -= dy * orbitSpeed;
            elevation = glm::clamp(elevation, 0.01f, float(M_PI) - 0.01f);
        }
        lastX = x;
        lastY = y;
        update();
    }
    void processMouseButton(int button, int action, int mods, GLFWwindow* win) {
        if (button == GLFW_MOUSE_BUTTON_LEFT || button == GLFW_MOUSE_BUTTON_MIDDLE) {
            if (action == GLFW_PRESS) {
                dragging = true;
                glfwGetCursorPos(win, &lastX, &lastY);
            } else if (action == GLFW_RELEASE) {
                dragging = false;
            }
        }
    }
    void processScroll(double xoffset, double yoffset) {
        radius -= yoffset * zoomSpeed;
        if (radius < 1.0f) radius = 1.0f;
        update();
    };
};
Camera camera;
vec3 sphericalToCartesian(float r, float theta, float phi){
        float x = r * sin(theta) * cos(phi);
        float y = r * cos(theta);
        float z = r * sin(theta) * sin(phi);
        return vec3(x, y, z);
    }
void generateParticles(int N) {
    particles.clear();
    for (int i = 0; i < N; ++i) {
        // --- get x, y, z, positions
        vec3 pos = sphericalToCartesian(
            sampleR(n, l, gen), 
            sampleTheta(l, m, gen), 
            samplePhi(n, l, m)
        );
        // --- color & add particle ---
        float r = length(pos);
        double theta = acos(pos.y / r);
        double phi = atan2(pos.z, pos.x);
        vec4 col = inferno(r, theta, phi, n, l, m) ;
        particles.emplace_back(pos, col);
    }
}

void buildQuantumCombos(int maxN = selectorMaxN) {
    quantumCombos.clear();
    for (int nValue = 1; nValue <= maxN; ++nValue) {
        for (int lValue = 0; lValue < nValue; ++lValue) {
            for (int mValue = -lValue; mValue <= lValue; ++mValue) {
                quantumCombos.push_back({nValue, lValue, mValue});
            }
        }
    }
}

int findComboIndex(int nValue, int lValue, int mValue) {
    for (size_t i = 0; i < quantumCombos.size(); ++i) {
        const QuantumCombo& combo = quantumCombos[i];
        if (combo.n == nValue && combo.l == lValue && combo.m == mValue) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

void syncCurrentComboIndex() {
    currentComboIndex = findComboIndex(n, l, m);
}

void stepComboSelection(int direction) {
    if (quantumCombos.empty()) return;

    if (currentComboIndex < 0) {
        syncCurrentComboIndex();
        if (currentComboIndex < 0) currentComboIndex = 0;
    }

    const int total = static_cast<int>(quantumCombos.size());
    currentComboIndex = (currentComboIndex + direction) % total;
    if (currentComboIndex < 0) currentComboIndex += total;

    const QuantumCombo& combo = quantumCombos[currentComboIndex];
    n = combo.n;
    l = combo.l;
    m = combo.m;
}

void updateWindowTitle(GLFWwindow* window) {
    if (!window) return;

    string title = currentElementLabel()
                 + " - n=" + to_string(n)
                 + " l=" + to_string(l)
                 + " m=" + to_string(m)
                 + " N=" + to_string(N);

    if (currentComboIndex >= 0 && !quantumCombos.empty()) {
        title += "  [" + to_string(currentComboIndex + 1)
              + "/" + to_string(quantumCombos.size()) + "]";
    }

    glfwSetWindowTitle(window, title.c_str());
}

using Glyph5x7 = array<uint8_t, 7>;

const Glyph5x7& glyphForChar(char c) {
    static const unordered_map<char, Glyph5x7> glyphs = {
        {' ', {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00}},
        {'-', {0x00, 0x00, 0x00, 0x1F, 0x00, 0x00, 0x00}},
        {'0', {0x0E, 0x11, 0x13, 0x15, 0x19, 0x11, 0x0E}},
        {'1', {0x04, 0x0C, 0x04, 0x04, 0x04, 0x04, 0x0E}},
        {'2', {0x0E, 0x11, 0x01, 0x02, 0x04, 0x08, 0x1F}},
        {'3', {0x1F, 0x01, 0x02, 0x06, 0x01, 0x11, 0x0E}},
        {'4', {0x02, 0x06, 0x0A, 0x12, 0x1F, 0x02, 0x02}},
        {'5', {0x1F, 0x10, 0x1E, 0x01, 0x01, 0x11, 0x0E}},
        {'6', {0x06, 0x08, 0x10, 0x1E, 0x11, 0x11, 0x0E}},
        {'7', {0x1F, 0x01, 0x02, 0x04, 0x08, 0x08, 0x08}},
        {'8', {0x0E, 0x11, 0x11, 0x0E, 0x11, 0x11, 0x0E}},
        {'9', {0x0E, 0x11, 0x11, 0x0F, 0x01, 0x02, 0x0C}},
        {'A', {0x0E, 0x11, 0x11, 0x1F, 0x11, 0x11, 0x11}},
        {'B', {0x1E, 0x11, 0x11, 0x1E, 0x11, 0x11, 0x1E}},
        {'C', {0x0E, 0x11, 0x10, 0x10, 0x10, 0x11, 0x0E}},
        {'D', {0x1C, 0x12, 0x11, 0x11, 0x11, 0x12, 0x1C}},
        {'E', {0x1F, 0x10, 0x10, 0x1E, 0x10, 0x10, 0x1F}},
        {'F', {0x1F, 0x10, 0x10, 0x1E, 0x10, 0x10, 0x10}},
        {'G', {0x0E, 0x11, 0x10, 0x17, 0x11, 0x11, 0x0F}},
        {'H', {0x11, 0x11, 0x11, 0x1F, 0x11, 0x11, 0x11}},
        {'I', {0x1F, 0x04, 0x04, 0x04, 0x04, 0x04, 0x1F}},
        {'J', {0x01, 0x01, 0x01, 0x01, 0x11, 0x11, 0x0E}},
        {'K', {0x11, 0x12, 0x14, 0x18, 0x14, 0x12, 0x11}},
        {'L', {0x10, 0x10, 0x10, 0x10, 0x10, 0x10, 0x1F}},
        {'M', {0x11, 0x1B, 0x15, 0x15, 0x11, 0x11, 0x11}},
        {'N', {0x11, 0x19, 0x15, 0x13, 0x11, 0x11, 0x11}},
        {'O', {0x0E, 0x11, 0x11, 0x11, 0x11, 0x11, 0x0E}},
        {'P', {0x1E, 0x11, 0x11, 0x1E, 0x10, 0x10, 0x10}},
        {'Q', {0x0E, 0x11, 0x11, 0x11, 0x15, 0x12, 0x0D}},
        {'R', {0x1E, 0x11, 0x11, 0x1E, 0x14, 0x12, 0x11}},
        {'S', {0x0F, 0x10, 0x10, 0x0E, 0x01, 0x01, 0x1E}},
        {'T', {0x1F, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04}},
        {'U', {0x11, 0x11, 0x11, 0x11, 0x11, 0x11, 0x0E}},
        {'V', {0x11, 0x11, 0x11, 0x11, 0x11, 0x0A, 0x04}},
        {'W', {0x11, 0x11, 0x11, 0x15, 0x15, 0x15, 0x0A}},
        {'X', {0x11, 0x11, 0x0A, 0x04, 0x0A, 0x11, 0x11}},
        {'Y', {0x11, 0x11, 0x0A, 0x04, 0x04, 0x04, 0x04}},
        {'Z', {0x1F, 0x01, 0x02, 0x04, 0x08, 0x10, 0x1F}}
    };

    static const Glyph5x7 fallback = {0x1F, 0x11, 0x02, 0x04, 0x04, 0x00, 0x04};
    const auto it = glyphs.find(c);
    return (it == glyphs.end()) ? fallback : it->second;
}

vector<float> buildTextVertices(const string& text, float originX, float originY, float pixelScale) {
    vector<float> vertices;
    vertices.reserve(text.size() * 7 * 5 * 12);

    float cursorX = originX;
    float cursorY = originY;
    const float charWidth = 5.0f * pixelScale;
    const float charHeight = 7.0f * pixelScale;
    const float spacing = pixelScale;

    for (char rawChar : text) {
        if (rawChar == '\n') {
            cursorX = originX;
            cursorY += charHeight + spacing * 2.0f;
            continue;
        }

        const char c = static_cast<char>(toupper(static_cast<unsigned char>(rawChar)));
        const Glyph5x7& glyph = glyphForChar(c);

        for (int row = 0; row < 7; ++row) {
            for (int col = 0; col < 5; ++col) {
                const bool bitOn = (glyph[row] >> (4 - col)) & 1;
                if (!bitOn) continue;

                const float x = cursorX + col * pixelScale;
                const float y = cursorY + row * pixelScale;
                const float x2 = x + pixelScale;
                const float y2 = y + pixelScale;

                vertices.insert(vertices.end(), {
                    x,  y,  x2, y,  x,  y2,
                    x2, y,  x2, y2, x,  y2
                });
            }
        }

        cursorX += charWidth + spacing;
    }

    return vertices;
}

struct Engine {
    GLFWwindow* window;
    int WIDTH = 800;
    int HEIGHT = 600;

    // renders vars
    GLuint sphereVAO, sphereVBO;
    int sphereVertexCount;
    GLuint shaderProgram;
    GLint modelLoc, viewLoc, projLoc, colorLoc;

    // text overlay vars
    GLuint textVAO = 0, textVBO = 0;
    GLuint textShaderProgram = 0;
    GLint textProjectionLoc = -1, textColorLoc = -1;

    // --- shaders ---
    const char* vertexShaderSource = R"glsl(
        #version 330 core
        layout(location=0) in vec3 aPos; uniform mat4 model; uniform mat4 view;
        uniform mat4 projection; out float lightIntensity;
        void main() { gl_Position = projection * view * model * vec4(aPos, 1.0);
            vec3 normal = normalize(aPos);
            vec3 lightDir = normalize(vec3(1.0, 1.0, 1.0));
            lightIntensity = max(dot(normal, lightDir), 0.5); // 0.2 is ambient light
        } )glsl";

    const char* fragmentShaderSource = R"glsl(
        #version 330 core
        in float lightIntensity; 
        out vec4 FragColor; 
        uniform vec4 objectColor;

        void main() { 
            // Increase the power to make the 'center-facing' spot tighter and brighter
            float glow = pow(lightIntensity, 2.0); 
            FragColor = vec4(objectColor.rgb , objectColor.a); 
        } )glsl";

    static bool checkShaderCompile(GLuint shader, const char* stage) {
        GLint ok = GL_FALSE;
        glGetShaderiv(shader, GL_COMPILE_STATUS, &ok);
        if (ok == GL_TRUE) return true;

        GLint logLen = 0;
        glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &logLen);
        string log(std::max(1, logLen), '\0');
        glGetShaderInfoLog(shader, logLen, nullptr, log.data());
        cerr << "Shader compile error (" << stage << "): " << log << '\n';
        return false;
    }

    static bool checkProgramLink(GLuint program) {
        GLint ok = GL_FALSE;
        glGetProgramiv(program, GL_LINK_STATUS, &ok);
        if (ok == GL_TRUE) return true;

        GLint logLen = 0;
        glGetProgramiv(program, GL_INFO_LOG_LENGTH, &logLen);
        string log(std::max(1, logLen), '\0');
        glGetProgramInfoLog(program, logLen, nullptr, log.data());
        cerr << "Program link error: " << log << '\n';
        return false;
    }

    void initTextRenderer() {
        const char* textVertexShader = R"glsl(
            #version 330 core
            layout(location = 0) in vec2 aPos;
            uniform mat4 projection;
            void main() {
                gl_Position = projection * vec4(aPos, 0.0, 1.0);
            }
        )glsl";

        const char* textFragmentShader = R"glsl(
            #version 330 core
            out vec4 FragColor;
            uniform vec4 textColor;
            void main() {
                FragColor = textColor;
            }
        )glsl";

        GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
        glShaderSource(vertexShader, 1, &textVertexShader, nullptr);
        glCompileShader(vertexShader);
        if (!checkShaderCompile(vertexShader, "text vertex")) {
            glDeleteShader(vertexShader);
            return;
        }

        GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
        glShaderSource(fragmentShader, 1, &textFragmentShader, nullptr);
        glCompileShader(fragmentShader);
        if (!checkShaderCompile(fragmentShader, "text fragment")) {
            glDeleteShader(vertexShader);
            glDeleteShader(fragmentShader);
            return;
        }

        textShaderProgram = glCreateProgram();
        glAttachShader(textShaderProgram, vertexShader);
        glAttachShader(textShaderProgram, fragmentShader);
        glLinkProgram(textShaderProgram);
        if (!checkProgramLink(textShaderProgram)) {
            glDeleteShader(vertexShader);
            glDeleteShader(fragmentShader);
            glDeleteProgram(textShaderProgram);
            textShaderProgram = 0;
            return;
        }

        textProjectionLoc = glGetUniformLocation(textShaderProgram, "projection");
        textColorLoc = glGetUniformLocation(textShaderProgram, "textColor");

        glDeleteShader(vertexShader);
        glDeleteShader(fragmentShader);

        glGenVertexArrays(1, &textVAO);
        glGenBuffers(1, &textVBO);
        glBindVertexArray(textVAO);
        glBindBuffer(GL_ARRAY_BUFFER, textVBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 12, nullptr, GL_DYNAMIC_DRAW);
        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);
        glBindVertexArray(0);
    }

    Engine() {
        if (!glfwInit()) {
            cerr << "GLFW init failed\n";
            exit(EXIT_FAILURE);
        }

#ifdef __APPLE__
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
        glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

        window = glfwCreateWindow(800, 600, "Atom Prob-Flow", NULL, NULL);
        if (!window) {
            cerr << "Failed to create GLFW window\n";
            glfwTerminate();
            exit(EXIT_FAILURE);
        }

        glfwMakeContextCurrent(window);
        glewExperimental = GL_TRUE;
        GLenum glewStatus = glewInit();
        if (glewStatus != GLEW_OK) {
            cerr << "GLEW init failed: " << reinterpret_cast<const char*>(glewGetErrorString(glewStatus)) << '\n';
            glfwDestroyWindow(window);
            glfwTerminate();
            exit(EXIT_FAILURE);
        }
        // Clear benign GL_INVALID_ENUM that can occur on core profiles after glewInit.
        glGetError();
        glEnable(GL_DEPTH_TEST);

        // Generate Sphere Vertices manually (like I did in the gravity sim)
        vector<float> vertices;
        float r = 0.05f; // Small sphere for particles
        int stacks = 10, sectors = 10;
        for(int i = 0; i <= stacks; ++i){
            float t1 = (float)i / stacks * M_PI;
            float t2 = (float)(i+1) / stacks * M_PI;
            for(int j = 0; j < sectors; ++j){
                float p1 = (float)j / sectors * 2 * M_PI;
                float p2 = (float)(j+1) / sectors * 2 * M_PI;
                auto getPos = [&](float t, float p) {
                    return vec3(r*sin(t)*cos(p), r*cos(t), r*sin(t)*sin(p));
                };
                vec3 v1 = getPos(t1, p1), v2 = getPos(t1, p2), v3 = getPos(t2, p1), v4 = getPos(t2, p2);
                vertices.insert(vertices.end(), {v1.x, v1.y, v1.z, v2.x, v2.y, v2.z, v3.x, v3.y, v3.z});
                vertices.insert(vertices.end(), {v2.x, v2.y, v2.z, v4.x, v4.y, v4.z, v3.x, v3.y, v3.z});
            }
        }
        sphereVertexCount = vertices.size() / 3;
        CreateVBOVAO(sphereVAO, sphereVBO, vertices);

        GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
        glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
        glCompileShader(vertexShader);
        if (!checkShaderCompile(vertexShader, "vertex")) {
            glfwDestroyWindow(window);
            glfwTerminate();
            exit(EXIT_FAILURE);
        }

        GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
        glShaderSource(fragmentShader, 1, &fragmentShaderSource, NULL);
        glCompileShader(fragmentShader);
        if (!checkShaderCompile(fragmentShader, "fragment")) {
            glfwDestroyWindow(window);
            glfwTerminate();
            exit(EXIT_FAILURE);
        }

        shaderProgram = glCreateProgram();
        glAttachShader(shaderProgram, vertexShader);
        glAttachShader(shaderProgram, fragmentShader);
        glLinkProgram(shaderProgram);
        if (!checkProgramLink(shaderProgram)) {
            glfwDestroyWindow(window);
            glfwTerminate();
            exit(EXIT_FAILURE);
        }

        // Get uniform locations
        modelLoc = glGetUniformLocation(shaderProgram, "model");
        viewLoc  = glGetUniformLocation(shaderProgram, "view");
        projLoc  = glGetUniformLocation(shaderProgram, "projection");
        colorLoc = glGetUniformLocation(shaderProgram, "objectColor");

        glDeleteShader(vertexShader);
        glDeleteShader(fragmentShader);
        initTextRenderer();
    }
    vec3 sphericalToCartesian(float r, float theta, float phi){
        float x = r * sin(theta) * cos(phi);
        float y = r * cos(theta);
        float z = r * sin(theta) * sin(phi);
        return vec3(x, y, z);
    }
    void CreateVBOVAO(GLuint& VAO, GLuint& VBO, const vector<float>& vertices) {
        glGenVertexArrays(1, &VAO);
        glGenBuffers(1, &VBO);
        glBindVertexArray(VAO);
        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_STATIC_DRAW);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);
    }
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
    void drawSpheres(vector<Particle>& particles) {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glUseProgram(shaderProgram); // Use our new shaded system

        int framebufferWidth = 0;
        int framebufferHeight = 0;
        glfwGetFramebufferSize(window, &framebufferWidth, &framebufferHeight);
        if (framebufferWidth <= 0 || framebufferHeight <= 0) return;

        glViewport(0, 0, framebufferWidth, framebufferHeight);
        const float aspect = static_cast<float>(framebufferWidth) / static_cast<float>(framebufferHeight);
        mat4 projection = perspective(radians(45.0f), aspect, 0.1f, 2000.0f);
        mat4 view = lookAt(camera.position(), camera.target, vec3(0, 1, 0)); 

        // Send view and projection to the shader
        glUniformMatrix4fv(viewLoc, 1, GL_FALSE, value_ptr(view));
        glUniformMatrix4fv(projLoc, 1, GL_FALSE, value_ptr(projection));

        glBindVertexArray(sphereVAO);

        for (auto& p : particles) {
            if (p.pos.x < 0 && p.pos.y > 0) continue;
            mat4 model = translate(mat4(1.0f), p.pos);
            model = scale(model, vec3(electron_r));
            glUniformMatrix4fv(modelLoc, 1, GL_FALSE, value_ptr(model));
            glUniform4f(colorLoc, p.color.r, p.color.g, p.color.b, p.color.a);
            
            glDrawArrays(GL_TRIANGLES, 0, sphereVertexCount);
        }
    }
    void drawTopLeftLabel(const string& text) {
        if (text.empty()) return;

        int framebufferWidth = 0;
        int framebufferHeight = 0;
        glfwGetFramebufferSize(window, &framebufferWidth, &framebufferHeight);
        if (framebufferWidth <= 0 || framebufferHeight <= 0) return;

        vector<float> textVertices = buildTextVertices(text, 16.0f, 16.0f, 3.0f);
        if (textVertices.empty()) return;

        vector<float> textVerticesClip;
        textVerticesClip.reserve((textVertices.size() / 2) * 3);
        for (size_t i = 0; i + 1 < textVertices.size(); i += 2) {
            const float px = textVertices[i];
            const float py = textVertices[i + 1];
            const float ndcX = (px / static_cast<float>(framebufferWidth)) * 2.0f - 1.0f;
            const float ndcY = 1.0f - (py / static_cast<float>(framebufferHeight)) * 2.0f;
            textVerticesClip.push_back(ndcX);
            textVerticesClip.push_back(ndcY);
            textVerticesClip.push_back(0.0f);
        }

        const bool depthWasEnabled = glIsEnabled(GL_DEPTH_TEST) == GL_TRUE;
        const bool blendWasEnabled = glIsEnabled(GL_BLEND) == GL_TRUE;

        glDisable(GL_DEPTH_TEST);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        glUseProgram(shaderProgram);
        const mat4 identity = mat4(1.0f);
        glUniformMatrix4fv(modelLoc, 1, GL_FALSE, value_ptr(identity));
        glUniformMatrix4fv(viewLoc, 1, GL_FALSE, value_ptr(identity));
        glUniformMatrix4fv(projLoc, 1, GL_FALSE, value_ptr(identity));
        glUniform4f(colorLoc, 1.0f, 1.0f, 1.0f, 1.0f);

        if (textVAO == 0 || textVBO == 0) {
            glGenVertexArrays(1, &textVAO);
            glGenBuffers(1, &textVBO);
        }

        glBindVertexArray(textVAO);
        glBindBuffer(GL_ARRAY_BUFFER, textVBO);
        glBufferData(GL_ARRAY_BUFFER, textVerticesClip.size() * sizeof(float), textVerticesClip.data(), GL_DYNAMIC_DRAW);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);
        glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(textVerticesClip.size() / 3));
        glBindVertexArray(0);

        if (!blendWasEnabled) glDisable(GL_BLEND);
        if (depthWasEnabled) glEnable(GL_DEPTH_TEST);
    }
    void setupCameraCallbacks() {
        glfwSetWindowUserPointer(window, &camera);
        glfwSetMouseButtonCallback(window, [](GLFWwindow* win, int button, int action, int mods) {
            ((Camera*)glfwGetWindowUserPointer(win))->processMouseButton(button, action, mods, win);
        });
        glfwSetCursorPosCallback(window, [](GLFWwindow* win, double x, double y) {
            ((Camera*)glfwGetWindowUserPointer(win))->processMouseMove(x, y);
        });
        glfwSetScrollCallback(window, [](GLFWwindow* win, double xoffset, double yoffset) {
            ((Camera*)glfwGetWindowUserPointer(win))->processScroll(xoffset, yoffset);
        });
        // Key callback: modify global quantum numbers / cycle valid combinations.
        glfwSetKeyCallback(window, [](GLFWwindow* win, int key, int scancode, int action, int mods) {
            if (!(action == GLFW_PRESS || action == GLFW_REPEAT)) return;

            bool shouldRegen = false;
            bool usedSelector = false;
            bool usedElementSelector = false;

            if (key == GLFW_KEY_UP) {
                stepElementSelection(1);
                shouldRegen = true;
                usedElementSelector = true;
            } else if (key == GLFW_KEY_DOWN) {
                stepElementSelection(-1);
                shouldRegen = true;
                usedElementSelector = true;
            } else if (key == GLFW_KEY_RIGHT) {
                stepComboSelection(1);
                shouldRegen = true;
                usedSelector = true;
            } else if (key == GLFW_KEY_LEFT) {
                stepComboSelection(-1);
                shouldRegen = true;
                usedSelector = true;
            } else if (key == GLFW_KEY_W) {
                n += 1;
                shouldRegen = true;
            } else if (key == GLFW_KEY_S) {
                n -= 1;
                if (n < 1) n = 1;
                shouldRegen = true;
            } else if (key == GLFW_KEY_E) {
                l += 1;
                shouldRegen = true;
            } else if (key == GLFW_KEY_D) {
                l -= 1;
                if (l < 0) l = 0;
                shouldRegen = true;
            } else if (key == GLFW_KEY_R) {
                m += 1;
                shouldRegen = true;
            } else if (key == GLFW_KEY_F) {
                m -= 1;
                shouldRegen = true;
            } else if (key == GLFW_KEY_T) {
                N += 100000;
                shouldRegen = true;
            } else if (key == GLFW_KEY_G) {
                N -= 100000;
                if (N < 1000) N = 1000;
                shouldRegen = true;
            }

            if (!shouldRegen) return;

            // Clamp to valid ranges
            if (l > n - 1) l = n - 1;
            if (l < 0) l = 0;
            if (m > l) m = l;
            if (m < -l) m = -l;

            electron_r = float(n) / 3.0f;
            syncCurrentComboIndex();
            if (!usedElementSelector) syncElementSelectionFromQuantum();
            cout << "Quantum numbers updated: n=" << n << " l=" << l << " m=" << m << " N=" << N << "\n";
            if (usedSelector && currentComboIndex >= 0) {
                cout << "Selector: " << (currentComboIndex + 1) << "/" << quantumCombos.size() << "\n";
            }
            if (usedElementSelector) {
                cout << "Element: " << currentElementLabel() << "\n";
            }
            generateParticles(N);
            updateWindowTitle(win);
        });
    }
};

struct Grid {
    Engine& engine;
    GLuint gridVAO, gridVBO;
    vector<float> vertices;
    Grid(Engine& engineRef) : engine(engineRef) {
        vertices = CreateGridVertices(500.0f, 2);
        engine.CreateVBOVAO(gridVAO, gridVBO, vertices.data(), vertices.size());
    }
    void Draw (GLint objectColorLoc) {
        glUseProgram(engine.shaderProgram);
        glUniform4f(objectColorLoc, 1.0f, 1.0f, 1.0f, 0.5f);
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

// ================= Main Loop ================= //
int main () {
    Engine engine;
    Grid grid(engine);
    buildQuantumCombos();
    buildElementMappings();
    setElementByIndex(0);
    syncCurrentComboIndex();

    GLint objectColorLoc = glGetUniformLocation(engine.shaderProgram, "objectColor");
    glUseProgram(engine.shaderProgram);
    engine.setupCameraCallbacks();

    // --- scale r for bigger orbitals ---
    electron_r = float(n) / 3.0f;

    // --- Sample particles ---
    N = 250000;
    generateParticles(N);
    updateWindowTitle(engine.window);
    cout << "Mapped known elements to (n,l,m): " << elementMappings.size() << " entries.\n";
    for (const ElementQuantum& element : elementMappings) {
        cout << setw(3) << element.atomicNumber << " "
             << element.name
             << " -> n=" << element.n
             << " l=" << element.l
             << " m=" << element.m << "\n";
    }
    cout << "Controls: UP/DOWN cycle elements, LEFT/RIGHT cycle valid (n,l,m) combinations.\n";

    float dt = 0.5f;
    cout << "Starting simulation..." << endl;
    while (!glfwWindowShouldClose(engine.window)) {
        grid.Draw(objectColorLoc);

        // ------ Update Probability current ------
        for (Particle& p : particles) {
            double r = length(p.pos);
            if (r > 1e-6) {
                double theta = acos(p.pos.y / r);
                p.vel = calculateProbabilityFlow(p, n, l, m);
                vec3 temp_pos = p.pos + p.vel * dt;
                double new_phi = atan2(temp_pos.z, temp_pos.x);
                p.pos = engine.sphericalToCartesian(r, theta, new_phi);
            }
        }
        // ------ Draw Particles ------
        engine.drawSpheres(particles);
        engine.drawTopLeftLabel(currentElementLabel());

        glfwSwapBuffers(engine.window);
        glfwPollEvents();
    }

    // --- close ---
    glfwDestroyWindow(engine.window);
    glfwTerminate();
    return 0;
}
