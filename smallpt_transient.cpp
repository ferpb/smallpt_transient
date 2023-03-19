// smallpt, a Path Tracer by Kevin Beason, 2009
// Adapted by Fernando Peña, 2023
// Comments based on presentation by David Cline

// Make : g++ -O3 -fopenmp smallpt_transient.cpp -o smallpt_transient
//        Remove "-fopenmp" for g++ version < 4.2
// Usage: time ./smallpt_transient 5000
//        convert -delay 5 -loop 0 [0-9]*.ppm image.gif

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

// Vector class, used for points, normals and colors
struct Vec
{
    double x, y, z; // position (x, y, z) or color (r,g,b)
    Vec(double x_ = 0, double y_ = 0, double z_ = 0)
    {
        x = x_;
        y = y_;
        z = z_;
    }
    Vec operator+(const Vec &b) const { return Vec(x + b.x, y + b.y, z + b.z); }
    Vec operator-(const Vec &b) const { return Vec(x - b.x, y - b.y, z - b.z); }
    Vec operator*(double b) const { return Vec(x * b, y * b, z * b); }
    Vec mult(const Vec &b) const { return Vec(x * b.x, y * b.y, z * b.z); }

    // Normalize the vector dividing by its length
    Vec &norm() { return *this = *this * (1 / sqrt(x * x + y * y + z * z)); }

    // Dot product of two vectors. It is equal to the cosine if both are normalized
    double dot(const Vec &b) const { return x * b.x + y * b.y + z * b.z; }

    // Cross products of two vectors
    Vec operator%(Vec &b) { return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x); } // cross
};

// Ray class
struct Ray
{
    // A ray is a parametric line with an origin (o) and a direction (d). Points along the ray
    // can be defined by varying the parameter (t): P(t) = o + t*d
    Vec o, d;
    Ray(Vec o_, Vec d_) : o(o_), d(d_) {}
};

// The surface reflection type (used in radiance())
enum Refl_t
{
    DIFF,
    SPEC,
    REFR
};

// smallpt only supports sphere objects!
struct Sphere
{
    // We can define a sphere by a center point (C) and a radius (r)
    // The implicit equation of the sphere in vector form is
    //     dot((P-C), (P-C)) - r^2 = 0
    double rad;  // radius
    Vec p, e, c; // position, emission, color
    Refl_t refl; // reflection type (DIFFuse, SPECular, REFRactive)
    Sphere(double rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_) : rad(rad_), p(p_), e(e_), c(c_), refl(refl_) {}

    // Ray-sphere intersection routine
    // Returns distance from ray origin to intersection or 0 if no hit
    double intersect(const Ray &r) const
    {
        // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
        Vec op = p - r.o;                            // p is the sphere center (C)
        double t, eps = 1e-4;                        // eps is a small fudge factor
        double b = op.dot(r.d);                      // 1/2 b from quadratic eq. setup
        double det = b * b - op.dot(op) + rad * rad; // (b^2-4ac)/4, a=1 because ray is normalized
        if (det < 0)                                 // ray missed the sphere
            return 0;
        else
            det = sqrt(det);
        return (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 0); // return smaller positive t
    }
};

// The hard coded scene
// Sphere(radius, position, emission, color, material)
Sphere spheres[] = {
    Sphere(1e5, Vec(1e5 + 1, 40.8, 81.6), Vec(), Vec(.75, .25, .25), DIFF),   // Left
    Sphere(1e5, Vec(-1e5 + 99, 40.8, 81.6), Vec(), Vec(.25, .25, .75), DIFF), // Rght
    Sphere(1e5, Vec(50, 40.8, 1e5), Vec(), Vec(.75, .75, .75), DIFF),         // Back
    Sphere(1e5, Vec(50, 40.8, -1e5 + 170), Vec(), Vec(), DIFF),               // Frnt
    Sphere(1e5, Vec(50, 1e5, 81.6), Vec(), Vec(.75, .75, .75), DIFF),         // Botm
    Sphere(1e5, Vec(50, -1e5 + 81.6, 81.6), Vec(), Vec(.75, .75, .75), DIFF), // Top
    Sphere(16.5, Vec(27, 16.5, 47), Vec(), Vec(1, 1, 1) * .999, SPEC),        // Mirr
    Sphere(16.5, Vec(73, 16.5, 78), Vec(), Vec(1, 1, 1) * .999, REFR),        // Glas
    Sphere(600, Vec(50, 681.6 - .27, 81.6), Vec(12, 12, 12), Vec(), DIFF)     // Light
};
int num_spheres = sizeof(spheres) / sizeof(Sphere);

// The output of the radiance function is a set of unbounded colors.
// We need to convert them to be between 0 and 255 to display them.

// Clamp function
inline double clamp(double x) { return x < 0 ? 0 : x > 1 ? 1
                                                         : x; }

// Converts float to integers to be saved in PPM file
// Applies a gamma correction of 2.2
inline int toInt(double x) { return int(pow(clamp(x), 1 / 2.2) * 255 + .5); }

// Routine to intersect rays with the scene of spheres
inline bool intersect(const Ray &r, double &t, int &id)
{
    // Check each sphere, one at a time. Keep the closest intersection
    double d;
    double inf = t = 1e20;

    for (int i = num_spheres; i--;)
        if ((d = spheres[i].intersect(r)) && d < t)
        {
            t = d;
            id = i;
        }
    return t < inf;
}

// Recursive routine that solves the rendering equation
// Computes the radiance estimate along ray r
//   r: the ray we are casting
//   depth: number of interactions in the current path
//   distance: distance travelled allong the current path
//   Xi: random number seed
//   E: wheter to include emissive color
//   return value: the vector with the radiance estimate
Vec radiance(const Ray &r, int &depth, double &distance, unsigned short *Xi, int E = 1)
{
    double t;   // distance to intersection
    int id = 0; // index of intersected object

    // Do intersection
    if (!intersect(r, t, id))
        return Vec(); // if miss, return black

    const Sphere &obj = spheres[id]; // the hit object

    distance = distance + t;

    // Surface properties
    Vec x = r.o + r.d * t;                // ray intersection point
    Vec n = (x - obj.p).norm();           // sphere normal
    Vec nl = n.dot(r.d) < 0 ? n : n * -1; // properly oriented surface normal
    Vec f = obj.c;                        // object color

    // Rusian Roulette
    // Stop the recursion randomly based on the surface reflectivity
    // using the maximum component (r,g,b) of the surface color.
    double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y
                                                        : f.z; // max refl

    // Don't do Russian Roulette before depth 5
    if (++depth > 5 || !p)
    {
        if (erand48(Xi) < p)
        {
            f = f * (1 / p);
        }
        else
        {
            return obj.e * E;
        }
    }

    if (obj.refl == DIFF)
    {
        // Ideal diffuse reflection

        // Sample a direction in the sphere
        double r1 = 2 * M_PI * erand48(Xi); // random angle around
        double r2 = erand48(Xi);
        double r2s = sqrt(r2); // random distance from center

        // Use normal to create orthonormal coordinate frame (w, y, v)
        Vec w = nl;                                                 // w = normal
        Vec u = ((fabs(w.x) > .1 ? Vec(0, 1) : Vec(1)) % w).norm(); // u is perpendicular to w
        Vec v = w % u;                                              // v is perpendicular to u and w

        // Sample unit hemisphere
        Vec d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm(); // d is random reflection ray

        return obj.e * E + f.mult(radiance(Ray(x, d), depth, distance, Xi, 1));
    }
    else if (obj.refl == SPEC)
    {
        // Ideal specular (mirror) reflection
        // Reflected ray: angle of incidence = angle of reflection
        return obj.e + f.mult(radiance(Ray(x, r.d - n * 2 * n.dot(r.d)), depth, distance, Xi));
    }

    // Otherwise, we have a dielectric (glass) surface

    // Ideal dielectric refraction

    // Glass is both reflective and refractive, we compute the reflected ray here
    Ray refl_ray(x, r.d - n * 2 * n.dot(r.d));

    // Determine if ray is entering of exiting glass
    bool into = n.dot(nl) > 0; // Ray from outside going in?

    double nc = 1;                         // Index of refraction for air
    double nt = 1.5;                       // Index of refraction for glass
    double nnt = into ? nc / nt : nt / nc; // nnt is 1/1.5 if ray goes air-glass
                                           // of 1.5, if goes glass-air
    double ddn = r.d.dot(nl);
    double cos2t;

    // Total internal reflection occurs when the light ray attempts to leave glass at a
    // too shallow angle
    // If the angle is to shallow (total internal reflection), all the light is reflected

    // If total internal reflection, reflect
    if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0) // Total internal reflection
        return obj.e + f.mult(radiance(refl_ray, depth, distance, Xi));

    // Otherwise, choose reflection or refraction using the Fresnel term

    // Compute the refracted ray
    Vec tdir = (r.d * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).norm();

    // R0 = reflectand at normal incidence based on index of refraction
    // c = 1 - cos(theta)
    // Re = fresnel reflectance
    double a = nt - nc;
    double b = nt + nc;
    double R0 = a * a / (b * b);
    double c = 1 - (into ? -ddn : tdir.dot(n));

    double Re = R0 + (1 - R0) * c * c * c * c * c;
    double Tr = 1 - Re;
    double P = .25 + .5 * Re;
    double RP = Re / P;
    double TP = Tr / (1 - P);

    // P = probability of reflecting
    // Russian roulette to sample reflection or refraction
    if (erand48(Xi) < P)
        return obj.e + f.mult(radiance(refl_ray, depth, distance, Xi) * RP);
    else
        return obj.e + f.mult(radiance(Ray(x, tdir), depth, distance, Xi) * TP);
}

// Main function, loops over image pixels, creates image, and saves it to a PPM file
int main(int argc, char *argv[])
{
    // Setup image and camera
    int w = 400, h = 400; // Image size

    // Transient rendering:
    // Temporal slices are defined using the optical distance (the distance that a ray travels across the scene)
    // instead of the time of flight to avoid extra divisions
    int start = 0, end = 1200, delta = 8;   // start and end distances, and step size betwen slices
    int num_slices = (end - start) / delta; // number of temporal slices
    bool continuous = true;                // Make lights emit a delta pulse or continous illumination

    // Read number of samples
    int samps = argc == 2 ? atoi(argv[1]) : 1;

    // Look from (pos) and gaze direction (dir)
    Ray cam(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).norm()); // camera (pos, dir)
    // Horizontal (right) camera direction (0.5135 defined the field of view angle)
    Vec cx = Vec(w * .5135 / h); // x direction increment (uses implicit 0 for y and z)
    // Vertical (up) vector of the camera
    Vec cy = (cx % cam.d).norm() * .5135; // y direction increment

    Vec r; // Used for colors of samples

    Vec *c = new Vec[w * h];                   // Image in stationary state
    Vec *bins = new Vec[num_slices * w * h];   // Images for each temporal slice
    int *counts = new int[num_slices * w * h]; // Number of samples that belong to each slice

    // Initialize counts
    for (int i = 0; i < num_slices * w * h; i++)
    {
        counts[i] = 0;
    }

// Run each iteration of the outer loop in its own thread
#pragma omp parallel for schedule(dynamic, 1) private(r) // OpenMP

    // Loop over all image pixels
    for (int y = 0; y < h; y++) // Loop over image rows
    {
        fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samps, 100. * y / (h - 1)); // print progress

        // Stores the state of the random number generator (erand48)
        // It's seeded using an arbitrary function of the row number to decorrelate
        // (at least visually) the sequences from row-to-row. In this way, the sequences
        // are deterministic and consistent from run to run, and independent of which thread
        // is executing and in what order the rows are executed.
        unsigned short Xi[3] = {0, 0, static_cast<unsigned short>(y * y * y)};

        for (int x = 0; x < w; x++) // Loop over image columns
        {
            Vec r = Vec();
            int i = (h - y - 1) * w + x;

            for (int s = 0; s < samps; s++)
            {
                // Anti-aliasing is done using supersampling inside each pixel, which removes all the jaggies
                // except arround the light.

                // r1 and r2 random values that determine the location of a sample within a pixel
                double r1 = 2 * erand48(Xi);
                double dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
                double r2 = 2 * erand48(Xi);
                double dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);

                // Ray direction and radiance
                // Compute ray direction using cam.d, cx and cy
                Vec d = cx * ((dx / 2 + x) / w - .5) +
                        cy * ((dy / 2 + y) / h - .5) + cam.d;

                // Use radiance function to estimate radiance
                int depth = 0;
                double distance = 0;
                Vec sample = radiance(Ray(cam.o + d * 140, d.norm()), depth, distance, Xi);
                r = r + sample * (1.0 / samps);

                if (distance < end)
                {
                    int t = distance / delta; // find time slice
                    int ti = t * h * w + i;
                    bins[ti] = bins[ti] + sample;
                    counts[ti] = counts[ti] + 1;

                    if (continuous)
                    {
                        // Save radiance in all future time slices
                        for (int tp = t + 1; tp < end / delta; tp++)
                        {
                            int ti = tp * h * w + i;
                            bins[ti] = bins[ti] + sample;
                            counts[ti] = counts[ti] + 1;
                        }
                    }
                }
            }

            // Add subpixel estimate
            c[i] = c[i] + Vec(clamp(r.x), clamp(r.y), clamp(r.z));
        }
    }

    fprintf(stderr, "\n");

    // Write stationary image in PPM format
    {
        FILE *f = fopen("image.ppm", "w");
        fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
        for (int i = 0; i < w * h; i++)
            fprintf(f, "%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
    }

    // Write transient images in PPM format
    {
        for (int t = 0; t < num_slices; t++)
        {
            char filename[10];
            sprintf(filename, "%04d.ppm", t);
            FILE *f = fopen(filename, "w");
            fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
            for (int i = 0; i < w * h; i++)
            {
                int ti = t * w * h + i;
                int n = counts[ti];
                Vec r;
                if (n > 0)
                    r = bins[ti] * (1.0 / n);
                fprintf(f, "%d %d %d ",
                        toInt(clamp(r.x)), toInt(clamp(r.y)), toInt(clamp(r.z)));
            }
        }
    }
}
