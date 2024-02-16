#include "rasterizer.h"

using namespace std;

namespace CGL {

    RasterizerImp::RasterizerImp(PixelSampleMethod psm, LevelSampleMethod lsm,
        size_t width, size_t height,
        unsigned int sample_rate) {
        this->psm = psm;
        this->lsm = lsm;
        this->width = width;
        this->height = height;
        this->sample_rate = sample_rate;

        sample_buffer.resize(width * height * sample_rate, Color::White);
    }

    // Used by rasterize_point and rasterize_line
    void RasterizerImp::fill_pixel(size_t x, size_t y, Color c) {
        // TODO: Task 2: You might need to this function to fix points and lines (such as the black rectangle border in test4.svg)
        // NOTE: You are not required to implement proper supersampling for points and lines
        // It is sufficient to use the same color for all supersamples of a pixel for points and lines (not triangles)



        double sample_size = sqrt(sample_rate);
        double sample_length = 1 / sample_size;
        for (double i = 0; i < 1; i += sample_length) {
            for (double j = 0; j < 1; j += sample_length) {
                
                sample_buffer[sample_rate * ((y * width) + x) + (j * sample_size) * sample_size + (sample_size * i)] = c;
            }
        }
    }

    // Rasterize a point: simple example to help you start familiarizing
    // yourself with the starter code.
    //
    void RasterizerImp::rasterize_point(float x, float y, Color color) {
        // fill in the nearest pixel
        int sx = (int)floor(x);
        int sy = (int)floor(y);

        // check bounds
        if (sx < 0 || sx >= width) return;
        if (sy < 0 || sy >= height) return;

        fill_pixel(sx, sy, color);
        return;
    }

    // Rasterize a line.
    void RasterizerImp::rasterize_line(float x0, float y0,
        float x1, float y1,
        Color color) {
        if (x0 > x1) {
            swap(x0, x1); swap(y0, y1);
        }

        float pt[] = { x0,y0 };
        float m = (y1 - y0) / (x1 - x0);
        float dpt[] = { 1,m };
        int steep = abs(m) > 1;
        if (steep) {
            dpt[0] = x1 == x0 ? 0 : 1 / abs(m);
            dpt[1] = x1 == x0 ? (y1 - y0) / abs(y1 - y0) : m / abs(m);
        }

        while (floor(pt[0]) <= floor(x1) && abs(pt[1] - y0) <= abs(y1 - y0)) {
            rasterize_point(pt[0], pt[1], color);
            pt[0] += dpt[0]; pt[1] += dpt[1];
        }
    }

    // Rasterize a triangle.
    void RasterizerImp::rasterize_triangle(float x0, float y0,
        float x1, float y1,
        float x2, float y2,
        Color color) {
        // TODO: Task 1: Implement basic triangle rasterization here, no supersampling

        // find dXi and dYi
        float dX0 = x1 - x0;
        float dY0 = y1 - y0;
        float dX1 = x2 - x1;
        float dY1 = y2 - y1;
        float dX2 = x0 - x2;
        float dY2 = y0 - y2;

        // resolve winding

        float wind = ((dX0 * dY1) - (dY0 * dX1) > 0) ? 1 : -1;

        // Super sampling counter
        double sample_size = sqrt(sample_rate);
        double sample_length = 1 / sample_size;

        // loop through points within bonding box
        
        for (double x = floor(min(x2, min(x0, x1))); x <= max(x2, max(x0, x1)); x++) {
            for (double y = floor(min(y2, min(y0, y1))); y <= max(y2, max(y0, y1)); y++) {
                for (double i = 0; i < 1; i += sample_length) {
                    for (double j = 0; j < 1; j += sample_length) {
                        double sample_x = x + i + (1 / (sample_size * 2));
                        double sample_y = y + j + (1 / (sample_size * 2));
                        if (wind * ((-(sample_x - x0) * dY0) + ((sample_y - y0) * dX0)) < 0) {
                            continue;
                        }
                        if (wind * ((-(sample_x - x1) * dY1) + ((sample_y - y1) * dX1)) < 0) {
                            continue;
                        }
                        if (wind * ((-(sample_x - x2) * dY2) + ((sample_y - y2) * dX2)) < 0) {
                            continue;
                        }
                        
                        sample_buffer[sample_rate * ((y * width) + x) + (j*sample_size) * sample_size + (sample_size*i)] = color;
                    }
                }
            }
        }



    }




    void RasterizerImp::rasterize_interpolated_color_triangle(float x0, float y0, Color c0,
        float x1, float y1, Color c1,
        float x2, float y2, Color c2)
    {
        // TODO: Task 4: Rasterize the triangle, calculating barycentric coordinates and using them to interpolate vertex colors across the triangle
        // Hint: You can reuse code from rasterize_triangle
        // find dXi and dYi
        float dX0 = x1 - x0;
        float dY0 = y1 - y0;
        float dX1 = x2 - x1;
        float dY1 = y2 - y1;
        float dX2 = x0 - x2;
        float dY2 = y0 - y2;

        // resolve winding

        float wind = ((dX0 * dY1) - (dY0 * dX1) > 0) ? 1 : -1;

        // Super sampling counter
        double sample_size = sqrt(sample_rate);
        double sample_length = 1 / sample_size;

        // loop through points within bonding box

        for (double x = floor(min(x2, min(x0, x1))); x <= max(x2, max(x0, x1)); x++) {
            for (double y = floor(min(y2, min(y0, y1))); y <= max(y2, max(y0, y1)); y++) {
                for (double i = 0; i < 1; i += sample_length) {
                    for (double j = 0; j < 1; j += sample_length) {
                        double sample_x = x + i + (1 / (sample_size * 2));
                        double sample_y = y + j + (1 / (sample_size * 2));
                        if (wind * ((-(sample_x - x0) * dY0) + ((sample_y - y0) * dX0)) < 0) {
                            continue;
                        }
                        if (wind * ((-(sample_x - x1) * dY1) + ((sample_y - y1) * dX1)) < 0) {
                            continue;
                        }
                        if (wind * ((-(sample_x - x2) * dY2) + ((sample_y - y2) * dX2)) < 0) {
                            continue;
                        }

                        // interpolation
                        // resolve winding
                        float xA = x0;
                        float yA = y0;
                        float xB = wind > 0 ? x1 : x2;
                        float yB = wind > 0 ? y1 : y2;
                        float xC = wind > 0 ? x2 : x1;
                        float yC = wind > 0 ? y2 : y1;

                        float a = (-(sample_x - xB) * (yC - yB) + (sample_y - yB) * (xC - xB)) /
                            (-(xA - xB) * (yC - yB) + (yA - yB) * (xC - xB));
                        float b = (-(sample_x - xC) * (yA - yC) + (sample_y - yC) * (xA - xC)) /
                            (-(xB - xC) * (yA - yC) + (yB - yC) * (xA - xC));
                        float r = 1 - a - b;
                        Color color = (a * c0 + b * c1 + r * c2);
                        sample_buffer[sample_rate * ((y * width) + x) + (j * sample_size) * sample_size + (sample_size * i)] = color;
                    }
                }
            }
        }


    }


    void RasterizerImp::rasterize_textured_triangle(float x0, float y0, float u0, float v0,
        float x1, float y1, float u1, float v1,
        float x2, float y2, float u2, float v2,
        Texture& tex)
    {
        // TODO: Task 5: Fill in the SampleParams struct and pass it to the tex.sample function.
        // TODO: Task 6: Set the correct barycentric differentials in the SampleParams struct.
        // Hint: You can reuse code from rasterize_triangle/rasterize_interpolated_color_triangle
        // find dXi and dYi
        float dX0 = x1 - x0;
        float dY0 = y1 - y0;
        float dX1 = x2 - x1;
        float dY1 = y2 - y1;
        float dX2 = x0 - x2;
        float dY2 = y0 - y2;

        // resolve winding

        float wind = ((dX0 * dY1) - (dY0 * dX1) > 0) ? 1 : -1;

        // Super sampling counter
        double sample_size = sqrt(sample_rate);
        double sample_length = 1 / sample_size;

        // loop through points within bonding box

        for (double x = floor(min(x2, min(x0, x1))); x <= max(x2, max(x0, x1)); x++) {
            for (double y = floor(min(y2, min(y0, y1))); y <= max(y2, max(y0, y1)); y++) {
                for (double i = 0; i < 1; i += sample_length) {
                    for (double j = 0; j < 1; j += sample_length) {
                        double sample_x = x + i + (1 / (sample_size * 2));
                        double sample_y = y + j + (1 / (sample_size * 2));
                        if (wind * ((-(sample_x - x0) * dY0) + ((sample_y - y0) * dX0)) < 0) {
                            continue;
                        }
                        if (wind * ((-(sample_x - x1) * dY1) + ((sample_y - y1) * dX1)) < 0) {
                            continue;
                        }
                        if (wind * ((-(sample_x - x2) * dY2) + ((sample_y - y2) * dX2)) < 0) {
                            continue;
                        }

                        // interpolation
                        // resolve winding
                        double xA = x0;
                        double yA = y0;
                        double xB = x1;
                        double yB = y1;
                        double xC = x2;
                        double yC = y2;
                        

                        // find barycentric coordinate
                        double a = (-(sample_x - xB) * (yC - yB) + (sample_y - yB) * (xC - xB)) /
                            (-(xA - xB) * (yC - yB) + (yA - yB) * (xC - xB));
                        double b = (-(sample_x - xC) * (yA - yC) + (sample_y - yC) * (xA - xC)) /
                            (-(xB - xC) * (yA - yC) + (yB - yC) * (xA - xC));
                        double r = 1 - a - b;

                        // get u v using the barycentric coordinate
                        double u, v;
                        u = a * u0 + b * u1 + r * u2;
                        v = a * v0 + b * v1 + r * v2;
                        
                        Vector2D uv(u, v);
                        Color color;
                        if (psm == P_NEAREST) {
                            color = tex.sample_nearest(uv);
                        }
                        else {
                            
                            color = tex.sample_bilinear(uv);
                        }
                        //std::cout << "uv: " << uv << std::endl;
                        sample_buffer[sample_rate * ((y * width) + x) + (j * sample_size) * sample_size + (sample_size * i)] = color;
                    }
                }
            }
        }



    }

    void RasterizerImp::set_sample_rate(unsigned int rate) {
        // TODO: Task 2: You may want to update this function for supersampling support


        this->sample_rate = rate;
        std::cout << "set sample rate: " << sample_rate << std::endl;

        this->sample_buffer.resize(width * height * sample_rate, Color::White);
    }


    void RasterizerImp::set_framebuffer_target(unsigned char* rgb_framebuffer,
        size_t width, size_t height)
    {
        // TODO: Task 2: You may want to update this function for supersampling support
          // no need to change

        this->width = width;
        this->height = height;
        this->rgb_framebuffer_target = rgb_framebuffer;


        this->sample_buffer.resize(this->width * this->height * this->sample_rate, Color::White);
    }


    void RasterizerImp::clear_buffers() {
        std::fill(rgb_framebuffer_target, rgb_framebuffer_target + 3 * width * height, 255);
        std::fill(sample_buffer.begin(), sample_buffer.end(), Color::White);
    }


    // This function is called at the end of rasterizing all elements of the
    // SVG file.  If you use a supersample buffer to rasterize SVG elements
    // for antialising, you could use this call to fill the target framebuffer
    // pixels from the supersample buffer data.
    //
    void RasterizerImp::resolve_to_framebuffer() {
        // TODO: Task 2: You will likely want to update this function for supersampling support
        // no need to change as I averaged the color in sample_buffer already


        double sample_size = sqrt(sample_rate);
        double sample_length = 1 / sample_size;



        for (int x = 0; x < width; ++x) {
            for (int y = 0; y < height; ++y) {
                Color col = Color();
                
                for (double i = 0; i < 1; i += sample_length) {
                    
                    for (double j = 0; j < 1; j += sample_length) {
                        /*std::cout << "debug: " << i << j << std::endl;*/
                        Color sample_col = sample_buffer[sample_rate * ((y * width) + x) + (j * sample_size) * sample_size + (sample_size * i)];
                        col.r += sample_col.r;
                        col.g += sample_col.g;
                        col.b += sample_col.b;
                    }
                    /*std::cout << "end of first inner loop: " << std::endl;*/
                }
                
                col *= (1.0 / (sample_rate));
                //std::cout << "color: " << col.b << std::endl;
                //col = Color::Black;
                for (int k = 0; k < 3; ++k) {
                    this->rgb_framebuffer_target[3 * (y * width + x) + k] = (&col.r)[k] * 255;
                }
            }
        }

    }

    Rasterizer::~Rasterizer() { }


}// CGL


