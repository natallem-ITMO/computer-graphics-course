#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <tuple>
#include <cassert>
#include <functional>
#include <complex>
#include <string>

static const double MY_PI = atan(1) * 4;

using uchar = unsigned char;

std::string input_file_name;
std::string output_file_name;
int result_width, result_height;
double dx, dy, input_gamma;
double B = 0;
double C = 0.5;
int input_scaling_way;

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

    bool read_header() {
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

bool program_error = false;
file_worker input(program_error);
file_worker output(program_error);

double to_linear_space(double d) {
    if (input_gamma == 0) {
        if (d <= 0.04045) {
            return d / 12.92;
        } else {
            return std::pow((d + 0.055) / 1.055, 2.4);
        }
    } else {
        return pow(d, input_gamma);
    }
}

double to_visual_space(double d) {
    if (input_gamma == 0) {
        if (d <= 0.0031308) {
            return d * 12.92;
        } else {
            return 1.055 * std::pow(d, 1 / 2.4) - 0.055;
        }
    } else {
        return pow(d, 1 / input_gamma);
    }
}

double nearest_neighbour_kernel(double d) {
    if (d > 0.5) return 0;
    if (d <= -0.5) return 0;
    return 1;
}

double bilinear_kernel(double d) {
    if (d >= 1 || d <= -1) return 0;
    if (d > 0) {
        return (1 - d);
    } else {
        return (1 + d);
    }
}

double lanczos3_kernel(double x) {
    if (x == 0) return 1;
    static const int a = 3;
    if (x >= -a && x < a)
        return (a * sin(MY_PI * x) * sin(MY_PI * x / a)) / (MY_PI * MY_PI * x * x);
    return 0;
}

double bc_spline_kernel(double x) {
    x = abs(x);
    double res = 0;
    if (x < 1) {
        res = (12 - 9 * B - 6 * C) * (x * x * x) + (-18 + 12 * B + 6 * C) * x * x + (6 - 2 * B);
    } else if (x >= 1 && x < 2) {
        res = (-B - 6 * C) * (x * x * x) + (6 * B + 30 * C) * (x * x) + (-12 * B - 48 * C) * x + (8 * B + 24 * C);
    }
    return res / 6.;
}

using modife_func_ptr_t = double (*)(double);

const std::vector<modife_func_ptr_t> kernel_functions = {
        &nearest_neighbour_kernel,
        &bilinear_kernel,
        &lanczos3_kernel,
        &bc_spline_kernel
};

const std::vector<double> radius_for_karnel = {
        1,
        1,
        3,
        2
};


using point = std::complex<double>;
using point_i = std::complex<int>;

double get_coordinate_by_pixel(int x, int length) {
    double center = (length / 2);
    double res = x - center;
    if (length % 2 == 0)
        res += 0.5;
    return res;
}

point get_coordinates_by_pixel(int x, int y, int width, int height) {
    return point(get_coordinate_by_pixel(x, width), get_coordinate_by_pixel(y, height));
}

point transform_to_input_coordinates(point c) {
    point delta(dx, -dy);
    c -= delta;
    return point(c.real() * (input.width / (double) result_width), c.imag() * (input.height / (double) result_height));
}

std::vector<double> find_nearest_coordinates(double x_input, int length_input, double coeff) {
    std::vector<double> result;
    double max, min, start, end;
    double radius = radius_for_karnel[input_scaling_way];
    if (coeff > 1) {
        radius *= coeff;
    }
    min = x_input - radius;
    max = x_input + radius;
    start = (min > 0) ? min + 1 : min - 1;
    end = (max > 0) ? max + 1 : max - 1;
    start = trunc(start);
    end = trunc(end);
    if (length_input % 2 == 0) {
        start -= 0.5;
        end += 0.5;
    }
    while (start <= end) {
        if (start <= max && start >= min) {
            result.push_back(start);
        }
        start += 1;
    }
    return result;
}

int get_pixel_by_coordinate(double x, int length) {
    int xx = x;
    xx = trunc(xx);
    if (x < 0) {
        if (x - trunc(xx) != 0)
            --xx;
    }
    xx += ((length - (length % 2)) / 2);
    return xx;
}

std::complex<int> get_pixel_by_coordinates(point x, int width, int height) {
    return std::complex<int>(get_pixel_by_coordinate(x.real(), width), get_pixel_by_coordinate(x.imag(), height));
}

std::vector<int> get_pixels_by_coordinate(const std::vector<double> &xs, int length) {
    std::vector<int> result;
    for (auto t : xs) {
        result.push_back(get_pixel_by_coordinate(t, length));
    }
    return result;
}

double get_color_with_kernel(const std::vector<double> &color_values, const std::vector<double> &distance_values,
                             double point) {
    double result = 0;
    double coeff_sum = 0;
    std::vector<double> percentage_values;
    for (int i = 0; i < color_values.size(); ++i) {
        double distance = distance_values[i] - point;
        double percentage = kernel_functions[input_scaling_way](distance);
        percentage_values.push_back(percentage);
        coeff_sum += percentage;
    }
    for (auto &p : percentage_values) {
        p /= coeff_sum;
    }
    for (int i = 0; i < color_values.size(); ++i) {
        double color = color_values[i];
        double percentage = percentage_values[i];
        result += color * percentage;
    }
    return result;
}

template<typename T>
T boarder(T x, T min, T max) {
    if (x < min) {
        return min;
    }
    if (x > max)
        return max;
    return x;
}

uchar get_pixel_color(int x, int y, int rgb = 0) {
    x = boarder<int>(x, 0, input.width - 1);
    y = boarder<int>(y, 0, input.height - 1);
    if (rgb) {
        return input.input_buffer[(x + y * input.width) * 3 + (rgb - 1)];
    }
    return input.input_buffer[x + y * input.width];
}

void transform_points(point &p, std::vector<double> &points, double coeff, bool is_width) {
    if (is_width) {
        p = point(p.real() / coeff, p.imag());
    } else {
        p = point(p.real(), p.imag() / coeff);
    }
    std::vector<double> new_points;
    for (auto t : points) {
        new_points.push_back(t / coeff);
    }
    points = new_points;
}

void get_color_for_result_pixel(int x, int y) {

    double scale_coeff_width = input.width / double(result_width);
    double scale_coeff_height = input.height / double(result_height);

    point result_coordinate = get_coordinates_by_pixel(x, y, result_width, result_height);
    point input_coordinate = transform_to_input_coordinates(result_coordinate);

    std::vector<double> neighbours_width = find_nearest_coordinates(input_coordinate.real(), input.width,
                                                                    scale_coeff_width);
    std::vector<double> neighbours_height = find_nearest_coordinates(input_coordinate.imag(), input.height,
                                                                     scale_coeff_height);
    std::vector<int> neighbours_pixels_width = get_pixels_by_coordinate(neighbours_width, input.width);
    std::vector<int> neighbours_pixels_height = get_pixels_by_coordinate(neighbours_height, input.height);
    if (scale_coeff_width > 1) {
        transform_points(input_coordinate, neighbours_width, scale_coeff_width, true);
    }
    if (scale_coeff_height > 1) {
        transform_points(input_coordinate, neighbours_height, scale_coeff_height, false);
    }
    int i_end = 3;
    if (input.is_P5) {
        i_end = 1;
    }
    for (int i = 1; i <= i_end; ++i) {
        std::vector<double> vertical_colors;
        for (auto y_pixel : neighbours_pixels_height) {
            std::vector<double> cur_colors;
            for (auto x_pixel : neighbours_pixels_width) {
                cur_colors.push_back(to_linear_space(get_pixel_color(x_pixel, y_pixel, (input.is_P5) ? 0 : i) / 255.));
            }
            double cur_color = get_color_with_kernel(cur_colors, neighbours_width, input_coordinate.real());
            vertical_colors.push_back(cur_color);
        }
        double final_color = get_color_with_kernel(vertical_colors, neighbours_height, input_coordinate.imag());
        final_color = boarder<double>(to_visual_space(final_color) * 255, 0, 255);
        output.output << (uchar) round(final_color);
    }
}

void process_scaling() {
    for (int h = 0; h < result_height; ++h) {
        for (int w = 0; w < result_width; ++w) {
            get_color_for_result_pixel(w, h);
        }
    }
}

bool read_arguments(int argc, char *argv[]) {
    static std::string usage = "Usage : lab3.exe <name_of_input_file> <name_of_output_files> <result_width> <result_height> <dx> <dy> <gamma> <scaling way> [<B> <C>]\n";
    if (argc != 9 && argc != 11) {
        std::cerr << "Incorrect arguments number\n" << usage;
        return false;
    }
    input_file_name = argv[1];
    output_file_name = argv[2];
    try {
        result_width = std::stoi(argv[3]);
        if (result_width < 1) {
            std::cerr << "Incorrect number for result width(value should be >= 1)\n" << usage;
            return false;
        }
        result_height = std::stoi(argv[4]);
        if (result_height < 1) {
            std::cerr << "Incorrect number for result height(value should be >= 1)\n" << usage;
            return false;
        }
        dx = std::stod(argv[5]);
        dy = std::stod(argv[6]);
        input_gamma = std::stod(argv[7]);
        input_scaling_way = std::stoi(argv[8]);
        if (input_scaling_way < 0 || input_scaling_way > 3) {
            std::cerr << "Incorrect number for scaling way(value should be in range [0..3])\n" << usage;
            return false;
        }
        if (argc == 11) {
            if (input_scaling_way != 3) {
                std::cerr << "Incorrect arguments number\n" << usage;
                return false;
            }
            B = std::stod(argv[9]);
            C = std::stod(argv[10]);
        }
    } catch (std::invalid_argument &ex) {
        std::cerr << "Cannot parse argument as number\n" << usage;
        return false;
    }
    return true;
}


int main(int argc, char *argv[]) {
    std::vector<std::pair<std::string, std::complex<int>>> names = {
            {"field.pgm",               std::complex<int>(1920, 1280)},
            {"color_field.ppm",         std::complex<int>(1280, 853)},
            {"test_gray.ppm",           std::complex<int>(258, 222)},
            {"silver_shining_test.ppm", std::complex<int>(128, 191)},
            {"color_changes_test.ppm",  std::complex<int>(1599, 1066)},
            {"gray.pgm",                std::complex<int>(400, 347)},
            {"res_gray.pgm",            std::complex<int>(1600, 1041)},
            {"faces.pgm",            std::complex<int>(500,500)},
    };
    bool testing = true;
    int i = 7;
    double coef_w = 3;
    double coef_h = 4;
    if (testing) {
        std::string name = names[i].first;
        int argct = 9;
        char **argvt = new char *[argct];
        argvt[0] = "lab3.exe";
        std::string name_file = "B:\\Projects\\GitProjects\\Graphics\\pictures\\input_pictures\\" + name;
        argvt[1] = const_cast<char *>(name_file.c_str());
        argvt[3] = const_cast<char *>((std::to_string(int(names[i].second.real() * coef_w))).c_str()); // result width
        argvt[4] = const_cast<char *>((std::to_string(int(names[i].second.imag() * coef_h))).c_str()); // result width
        argvt[5] = "0"; // dx
        argvt[6] = "0"; // dy
        argvt[7] = "0"; // gamma
        argvt[8] = "0"; // scaling way
        /*
 <способ_масштабирования>:
0 – ближайшая точка (метод ближайшего соседа);
1 – билинейное;
2 – Lanczos3;
3 – BC-сплайны. Для этого способа могут быть указаны ещё два параметра: B и C, по умолчанию 0 и 0.5 (Catmull-Rom).
         */
        std::string name_file_out =
                "B:\\Projects\\GitProjects\\Graphics\\pictures\\output_pictures\\lab_4_ex\\" +
                std::string(argvt[3]) + "_" +
                std::string(argvt[4]) + "_" +
                std::string(argvt[5]) + "_" +
                std::string(argvt[6]) + "_" +
                std::string(argvt[7]) + "_" +
                std::string(argvt[8]) + "_" +
                name;
        argvt[2] = const_cast<char *>(name_file_out.c_str());
        if (!read_arguments(argct, argvt)) {
            return 1;
        }
    } else {
        if (!read_arguments(argc, argv)) {
            return 1;
        }
    }

    input.input_file_name = input_file_name;
    if (!input.open_input_file()) {
        return 1;
    }
    if (!input.read_header()) {
        return 1;
    }
    if (!input.read_buffer()) {
        return 1;
    }
    output.output_file_name = output_file_name;
    if (!output.open_output_file()) {
        return 1;
    }
    output.set_header_for_output_file(input.is_P5, result_width, result_height, 255);
    output.write_header();
    process_scaling();
    return 0;
}
