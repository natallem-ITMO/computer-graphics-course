#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <tuple>
#include <set>
#include <map>
#include <iomanip>
#include <cassert>

const std::vector<std::string> color_spaces = {"RGB", "HSL", "HSV", "YCbCr.601", "YCbCr.709", "YCoCg", "CMY"};
struct pixel_d;
using modife_func_ptr_t = pixel_d(*)(pixel_d);

int input_color_space = -1;
int output_color_space = -1;
int number_of_input_files = 0;
std::string input_file_name;
int number_of_output_files = 0;
std::string output_file_name;

struct pixel_ch {
    unsigned char x, y, z;

    pixel_ch() = default;

    pixel_ch(unsigned char x, unsigned char y, unsigned char z) : x(x), y(y), z(z) {}
};

struct pixel_d {
    double x, y, z;

    pixel_d(pixel_ch p) {
        x = double(p.x) / 255.;
        y = double(p.y) / 255.;
        z = double(p.z) / 255.;
    }

    pixel_d(double x, double y, double z) : x(x), y(y), z(z) {}

    static double constraint(double x) {
        if (x < 0) return 0;
        if (x > 1) return 1;
        return x;
    }

    operator pixel_ch() {
        unsigned char x_ = round(constraint(x) * 255.);
        unsigned char y_ = round(constraint(y) * 255.);
        unsigned char z_ = round(constraint(z) * 255.);
        return pixel_ch{x_, y_, z_};
    }
};

pixel_d from_RGB_to_RGB(pixel_d rgb) {
    return rgb;
}

pixel_d from_YCbCr_to_RGB(pixel_d ycbcr, bool is_709 = false) {
    double Kr = 0.299;
    double Kg = 0.587;
    if (is_709) {
        Kr = 0.2126;
        Kg = 0.7152;
    }
    double Kb = 1 - Kr - Kg;
    double Y = ycbcr.x;
    double Cb = ycbcr.y - 0.5;
    double Cr = ycbcr.z - 0.5;
    double R = Y + (2. - 2. * Kr) * Cr;
    double G = Y - Kb / Kg * (2. - 2. * Kb) * Cb - Kr / Kg * (2. - 2. * Kr) * Cr;
    double B = Y + (2. - 2. * Kb) * Cb;
    return pixel_d(R, G, B);
}

pixel_d from_RGB_to_YCbCr(pixel_d rgb, bool is_709 = false) {
    double Kr = 0.299;
    double Kg = 0.587;
    if (is_709) {
         Kr = 0.2126;
         Kg = 0.7152;
    }
    double Kb = 1 - Kr - Kg;
    double R = rgb.x;
    double G = rgb.y;
    double B = rgb.z;
    double Y = Kr * R + Kg * G + Kb * B;
    double Cb = -0.5 * Kr / (1 - Kb) * R - 0.5 * Kg / (1 - Kb) * G + 0.5 * B;
    double Cr = 0.5 * R - 0.5 * Kg / (1 - Kr) * G - 0.5 * Kb / (1 - Kr) * B;
    return pixel_d(Y, Cb + 0.5, Cr + 0.5);
}

pixel_d from_RGB_to_YCbCr_601(pixel_d rgb) {
    return from_RGB_to_YCbCr(rgb);
}

pixel_d from_RGB_to_YCbCr_709(pixel_d rgb) {
    return from_RGB_to_YCbCr(rgb, true);
}

pixel_d from_YCbCr_601_to_RGB(pixel_d rgb) {
    return from_YCbCr_to_RGB(rgb, false);
}

pixel_d from_YCbCr_709_to_RGB(pixel_d rgb) {
    return from_YCbCr_to_RGB(rgb, true);
}

pixel_d from_YCoCg_to_RGB(pixel_d ycocg) {
    double Y = ycocg.x;
    double Co = ycocg.y - 0.5;
    double Cg = ycocg.z - 0.5;
    double R = Y + Co - Cg;
    double G = Y + Cg;
    double B = Y - Co - Cg;
    return pixel_d(R, G, B);
}

pixel_d from_RGB_to_YCoCg(pixel_d rgb) {
    double R = rgb.x;
    double G = rgb.y;
    double B = rgb.z;
    double Y = 0.25 * R + 0.5 * G + 0.25 * B;
    double Co = 0.5 * R - 0.5 * B;
    double Cg = -0.25 * R + 0.5 * G - 0.25 * B;
    return pixel_d(Y, Co + 0.5, Cg + 0.5);
}

pixel_d to_CMY(pixel_d p) {
    return pixel_d(1 - p.x, 1 - p.y, 1 - p.z);
}

auto mod = [](double a, double b) {
    while (a >= b) {
        a -= b;
    }
    return a;
};

pixel_d from_HSL_to_RGB(pixel_d hsl) {
    double H = hsl.x * 360;
    double S = hsl.y;
    double L = hsl.z;
    double H_ = H / 60;
    double C = (1 - std::abs(2 * L - 1)) * S;
    double X = C * (1 - std::abs(mod(H_, 2.) - 1));
    double m = L - C / 2;
    std::vector<double> rgb;
    if (H_ >= 0 && H_ < 1) {
        rgb = {C, X, 0};
    }
    if (H_ >= 1 && H_ < 2) {
        rgb = {X, C, 0};
    }
    if (H_ >= 2 && H_ < 3) {
        rgb = {0, C, X};
    }
    if (H_ >= 3 && H_ < 4) {
        rgb = {0, X, C};
    }
    if (H_ >= 4 && H_ < 5) {
        rgb = {X, 0, C};
    }
    if (H_ >= 5) {
        rgb = {C, 0, X};
    }
    return pixel_d(rgb[0] + m, rgb[1] + m, rgb[2] + m);
}

pixel_d from_RGB_to_HSL(pixel_d rgb) {
    double R = rgb.x;
    double G = rgb.y;
    double B = rgb.z;
    double V;
    double X_max = V = std::max(std::max(R, G), B);
    double X_min = std::min(std::min(R, G), B);
    double C = X_max - X_min;
    double L = (X_max + X_min) / 2;
    double H;
    if (std::abs(C) < 1e-9) {
        H = 0;
    } else if (V == R) {
        H = 60. * (mod(((G - B) / C) + 6, 6.));
    } else if (V == G) {
        H = 60. * (2 + (B - R) / C);
    } else if (V == B) {
        H = 60. * (4 + (R - G) / C);
    }
    double S = (C == 0) ? 0 : C / (1 - std::abs(2 * L - 1));
    return pixel_d(H / 360., S, L);
}

pixel_d from_HSV_to_RGB(pixel_d hsl) {
    double H = hsl.x * 360;
    double S = hsl.y;
    double V = hsl.z;
    double H_ = H / 60;
    double C = V * S;
    double X = C * (1 - std::abs(mod(H_, 2.) - 1));

    std::vector<double> rgb;
    if (H_ >= 0 && H_ < 1) {
        rgb = {C, X, 0};
    }
    if (H_ >= 1 && H_ < 2) {
        rgb = {X, C, 0};
    }
    if (H_ >= 2 && H_ < 3) {
        rgb = {0, C, X};
    }
    if (H_ >= 3 && H_ < 4) {
        rgb = {0, X, C};
    }
    if (H_ >= 4 && H_ < 5) {
        rgb = {X, 0, C};
    }
    if (H_ >= 5) {
        rgb = {C, 0, X};
    }
    double m = V - C;
    return pixel_d(rgb[0] + m, rgb[1] + m, rgb[2] + m);
}

pixel_d from_RGB_to_HSV(pixel_d rgb) {
    double R = rgb.x;
    double G = rgb.y;
    double B = rgb.z;
    double V;
    double X_max = V = std::max(std::max(R, G), B);
    double X_min = std::min(std::min(R, G), B);
    double C = X_max - X_min;
    double L = (X_max + X_min) / 2;
    double H;
    if (std::abs(C) < 1e-9) {
        H = 0;
    } else if (V == R) {
        H = 60. * (mod(((G - B) / C) + 6, 6.));
    } else if (V == G) {
        H = 60. * (2 + (B - R) / C);
    } else if (V == B) {
        H = 60. * (4 + (R - G) / C);
    }
    double S = V ? C / V : 0;
    return pixel_d(H / 360., S, V);
}

const std::vector<modife_func_ptr_t> modification_to_rgb = {
        &from_RGB_to_RGB,
        &from_HSL_to_RGB,
        &from_HSV_to_RGB,
        &from_YCbCr_601_to_RGB,
        &from_YCbCr_709_to_RGB,
        &from_YCoCg_to_RGB,
        &to_CMY,
};
const std::vector<modife_func_ptr_t> modification_from_rgb = {
        &from_RGB_to_RGB,
        &from_RGB_to_HSL,
        &from_RGB_to_HSV,
        &from_RGB_to_YCbCr_601,
        &from_RGB_to_YCbCr_709,
        &from_RGB_to_YCoCg,
        &to_CMY,
};

struct file_worker {
    file_worker(bool &error) : error(error) {}

    file_worker(const std::string &inputFileName, const std::string &outputFileName, bool &error) : input_file_name(
            inputFileName), output_file_name(outputFileName), error(error) {};

    bool open_input_file() {
        assert(!input_file_name.empty());
        input.open(input_file_name, std::ios::binary);
        if (!input.is_open()) {
            std::cerr << "Failed to open input file \"" << input_file_name << "\"\n";
            error = true;
            return false;
        }
        return true;
    }

    bool open_output_file() {
        output.open(output_file_name, std::ios::binary);
        if (!output.is_open()) {
            std::cerr << "Failed to open output file \"" << output_file_name << "\"\n";
            error = true;
            return false;
        }
        return true;
    }

    void set_header_for_output_file(bool is_P5_, size_t w, size_t h, size_t max_color) {
        is_P5 = is_P5_;
        width = w;
        height = h;
        max_color_size = max_color;
    }

    ~file_worker() {
        delete[] input_buffer;
        if (input.is_open()) {
            input.close();
        }
        if (output.is_open()) {
            output.close();
        }
    }

    bool read_header(bool need_to_be_P5 = true) {
        std::string format;
        input >> format;
        if (format != "P5" && format != "P6") {
            std::cerr << "Input image in file \"" << input_file_name << "\" is not in format P5 or P6\n";
            error = true;
            return false;
        }
        if (format != "P5") {
            is_P5 = false;
        }
        if (!is_P5 && need_to_be_P5) {
            std::cerr << "Input image in file \"" << input_file_name << "\" should be in pgm format\n";
            error = true;
            return false;
        }
        if (is_P5 && !need_to_be_P5) {
            std::cerr << "Input image in file \"" << input_file_name << "\" should be in ppm format\n";
            error = true;
            return false;
        }

        input >> width;
        if (width == 0) {
            std::cerr << "Width of input image in file \"" << input_file_name << "\" is incorrect (= 0)\n";
            error = true;
            return false;
        }
        input >> height;
        if (height == 0) {
            std::cerr << "Height of input image in file \"" << input_file_name << "\" is incorrect (= 0)\n";
            error = true;
            return false;
        }
        int temp;
        input >> temp;
        max_color_size = temp;
        if (temp <= 0) {
            std::cerr << "Maximum color number of input image in file \"" << input_file_name
                      << "\" is incorrect (<= 0)\n";
            error = true;
            return false;
        }
        if (max_color_size != 255) {
            std::cerr << "Maximum color number of input image in file \"" << input_file_name
                      << "\" should be 255\n";
            error = true;
            return false;
        }
        input.get();
        size = width * height;
        number_of_pixels = (is_P5) ? size : size * 3;
        return true;
    }

    bool compare_with_other_input_file(const file_worker &other) const {
        if ((other.width == width) && (other.height == height) && (other.max_color_size == max_color_size) &&
            (other.is_P5 == is_P5)) {
            return true;
        } else {
            error = true;
            return false;
        }
    }

    bool read_buffer() {
        try {
            input_buffer = new char[number_of_pixels];
            input.read(input_buffer, number_of_pixels);
            if (input.gcount() != number_of_pixels) {
                size_t d = input.gcount();
                std::cerr << "Couldn't read all bytes of input image from file \'" << input_file_name
                          << "\' to buffer\n";
                error = true;
                return false;
            }
        } catch (const std::bad_alloc &e) {
            std::cerr << "Allocation for input image from file \'" << input_file_name << "\' buffer failed: "
                      << e.what() << '\n';
            error = true;
            return false;
        }
        return true;
    }

    void write_header() {
        assert(output.is_open());
        output << (is_P5 ? "P5" : "P6") << '\n';
        output << std::to_string(width) + " " + std::to_string(height) << '\n';
        output << std::to_string(max_color_size) << '\n';
    }

    std::string input_file_name;
    std::string output_file_name;

    bool &error;

    char *input_buffer = nullptr;
    std::ifstream input;
    std::ofstream output;

    size_t width;
    size_t height;
    size_t max_color_size;
    size_t size;
    size_t number_of_pixels;

    bool is_P5 = true;
};

bool error = false;
file_worker input_1(error);
file_worker input_2(error);
file_worker input_3(error);
file_worker output_1(error);
file_worker output_2(error);
file_worker output_3(error);

bool read_arguments(int argc, char *argv[]) {
    static std::string usage = "Usage : lab2.exe -f <input_color_space> -t <output_color_space> "
                               "-i <number_of_input_files> <name_of_input_file> "
                               "-o <number_of_output_files> <name_of_output_files>\n";
    if (argc != 11) {
        std::cerr << "Incorrect arguments number\n" << usage;
        return false;
    }
    int index = 1;
    while (index != 11) {
        std::string mode = argv[index++];
        if (index == 11) {
            std::cerr << "Incorrect arguments\n" << usage;
            return false;
        }
        if (mode == "-f") {
            std::string input_space = argv[index++];
            for (int i = 0; i < color_spaces.size(); ++i) {
                if (color_spaces[i] == input_space) {
                    input_color_space = i;
                    break;
                }
            }
            if (input_color_space == -1) {
                std::cerr << "Incorrect color space name \'" << input_space << "\'\n" << usage;
                return false;
            }
        } else if (mode == "-t") {
            std::string output_space = argv[index++];
            for (int i = 0; i < color_spaces.size(); ++i) {
                if (color_spaces[i] == output_space) {
                    output_color_space = i;
                    break;
                }
            }
            if (output_color_space == -1) {
                std::cerr << "Incorrect color space name \'" << output_space << "\'\n" << usage;
                return false;
            }
        } else if (mode == "-i") {
            try {
                number_of_input_files = std::stoi(argv[index++]);
                if (number_of_input_files != 1 && number_of_input_files != 3) {
                    std::cerr << "Incorrect number of input files\n" << usage;
                    return false;
                }
            } catch (std::invalid_argument &ex) {
                std::cerr << "Cannot parse argument as integer\n" << usage;
                return false;
            }
            if (index == 11) {
                std::cerr << "Incorrect arguments\n" << usage;
                return false;
            }
            input_file_name = argv[index++];
        } else if (mode == "-o") {
            try {
                number_of_output_files = std::stoi(argv[index++]);
                if (number_of_output_files != 1 && number_of_output_files != 3) {
                    std::cerr << "Incorrect number of output files\n" << usage;
                    return false;
                }
            } catch (std::invalid_argument &ex) {
                std::cerr << "Cannot parse argument as integer\n" << usage;
                return false;
            }
            if (index == 11) {
                std::cerr << "Incorrect arguments\n" << usage;
                return false;
            }
            output_file_name = argv[index++];
        } else {
            std::cerr << "Incorrect arguments \'" << mode << "\'\n" << usage;
            return false;
        }
    }
    if (input_file_name.empty() || output_file_name.empty() || input_color_space == -1 || output_color_space == -1 ||
        !number_of_input_files || !number_of_output_files) {
        std::cerr << "Lack of arguments\n" << usage;
        return false;
    }
    return true;
}

std::pair<std::string, std::string> parse_file_name(std::string &str) {
    int index_before_last_point = -2;
    for (int i = str.size() - 1; i >= 0; i--) {
        if (str[i] == '.') {
            index_before_last_point = i - 1;
            break;
        }
    }
    if (index_before_last_point != -2) {
        if (index_before_last_point == -1) {
            return std::make_pair("", str);
        } else {
            return std::make_pair(str.substr(0, index_before_last_point + 1),
                                  str.substr(index_before_last_point + 1, str.size()));
        }
    } else {
        std::cerr << "Incorrect file name \'" << str << "\'\n";
        return {"", ""};
    }
}

pixel_ch read_pixel(size_t i) {
    if (number_of_input_files == 3) {
        return pixel_ch((unsigned char) input_1.input_buffer[i], (unsigned char) input_2.input_buffer[i],
                        (unsigned char) input_3.input_buffer[i]);
    } else {
        return pixel_ch((unsigned char) input_1.input_buffer[3 * i], (unsigned char) input_1.input_buffer[3 * i + 1],
                        (unsigned char) input_1.input_buffer[3 * i + 2]);
    }
}

void write_pixel(pixel_ch ch) {
    if (number_of_output_files == 3) {
        output_1.output << ch.x;
        output_2.output << ch.y;
        output_3.output << ch.z;
    } else {
        output_1.output << ch.x << ch.y << ch.z;
    }
}

void transform() {
    for (size_t i = 0; i < input_1.size; i++) {
        pixel_ch input = read_pixel(i);
        pixel_d in_rgb = modification_to_rgb[input_color_space](input);
        pixel_d output = modification_from_rgb[output_color_space](in_rgb);
        write_pixel(output);
    }
}

int main(int argc, char *argv[]) {
    if (!read_arguments(argc, argv)) {
        return 1;
    }

    if (number_of_input_files == 3) {
        auto p_input = parse_file_name(input_file_name);
        if (p_input.second.empty()) {
            return 1;
        }
        input_1.input_file_name = p_input.first + "_1" + p_input.second;
        input_2.input_file_name = p_input.first + "_2" + p_input.second;
        input_3.input_file_name = p_input.first + "_3" + p_input.second;
        if (!input_1.open_input_file() || !input_2.open_input_file() || !input_3.open_input_file()) {
            return 1;
        }
        if (!input_1.read_header() || !input_2.read_header() || !input_3.read_header()) {
            return 1;
        }
        if (!input_1.compare_with_other_input_file(input_2) || !input_1.compare_with_other_input_file(input_3)) {
            return 1;
        }
        if (!input_1.read_buffer() || !input_2.read_buffer() || !input_3.read_buffer()) {
            return 1;
        }
    } else {
        input_1.input_file_name = input_file_name;
        if (!input_1.open_input_file() || !input_1.read_header(false) || !input_1.read_buffer()) {
            if (error) {
                return 1;
            }
        }
    }

    if (number_of_output_files == 3) {
        auto p_output = parse_file_name(output_file_name);
        if (p_output.second.empty()) {
            return 1;
        }
        output_1.output_file_name = p_output.first + "_1" + p_output.second;
        output_2.output_file_name = p_output.first + "_2" + p_output.second;
        output_3.output_file_name = p_output.first + "_3" + p_output.second;
        if (!output_1.open_output_file() || !output_2.open_output_file() || !output_3.open_output_file()) {
            return 1;
        }
        output_1.set_header_for_output_file(true, input_1.width, input_1.height, 255);
        output_2.set_header_for_output_file(true, input_1.width, input_1.height, 255);
        output_3.set_header_for_output_file(true, input_1.width, input_1.height, 255);

        output_1.write_header();
        output_2.write_header();
        output_3.write_header();
    } else {
        output_1.output_file_name = output_file_name;
        if (!output_1.open_output_file()) {
            return 1;
        }
        output_1.set_header_for_output_file(false, input_1.width, input_1.height, 255);
        output_1.write_header();
    }
    transform();

    return 0;
}
