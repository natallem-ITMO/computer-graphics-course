#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <tuple>
#include <set>
#include <map>
#include <iomanip>
#include <cassert>
#include <numeric>
#include <random>
#include <functional>

using uchar = unsigned char;
int num_of_colors;
std::vector<uchar> num_to_byte;
std::string input_file_name;
std::string output_file_name;
int input_gradient;
int input_dithering;
int input_bits;
double input_gamma;

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

template<typename T>
struct buffer {

    buffer() = default;

    bool alloc(int ww, int hh) {
        w = ww;
        h = hh;
        try {
            buf = new T[ww * hh];
        } catch (const std::bad_alloc &e) {
            std::cerr << "Cannot alloc buffer for computation\n";
            return false;
        }
        return true;
    }

    T &get_ref(int x, int y) {
        return buf[x + y * w];
    }

    void fill_gradient() {
        for (int x = 0; x < w; ++x) {
            double val = (num_of_colors - 1) * double(x) / (w - 1);
            for (int y = 0; y < h; ++y) {
                get_ref(x, y) = val;
            }
        }
    }

    void fill_with_data() {
        for (int y = 0; y < h; ++y) {
            for (int x = 0; x < w; ++x) {
                get_ref(x, y) = ((uchar) input.input_buffer[x + y * w]) / 255. * (num_of_colors - 1);
            }
        }
    }

    void fill_zeros() {
        for (int x = 0; x < w; ++x) {
            for (int y = 0; y < h; ++y) {
                get_ref(x, y) = 0;
            }
        }
    };

    ~buffer() {
        delete buf;
    }

    T *buf = nullptr;
    int w = 0, h = 0;
};

buffer<double> gradient_buffer;
buffer<double> error_buffer;

double to_liniar_space(double d) {
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

void no_dithering() {
    for (int i = 0; i < input.size; i++) {
        output.output << (uchar) num_to_byte[(int) round(gradient_buffer.buf[i])];
    }
}

std::pair<uchar, uchar> give_neighbours(double &ref) {
    int prev = trunc(ref);
    if (prev == (num_of_colors - 1)) {
        return {num_to_byte[num_of_colors - 2], num_to_byte[num_of_colors - 1]};
    }
    return {num_to_byte[prev], num_to_byte[prev + 1]};
}

void dithering_error_not_saving(std::function<double(int, int)> error_get_func) {
    for (int y = 0; y < input.height; y++) {
        for (int x = 0; x < input.width; ++x) {
            double cur_val = gradient_buffer.get_ref(x, y);
            auto[prev, next] = give_neighbours(cur_val);
            double prev_gamma = to_liniar_space(prev / 255.);
            double next_gamma = to_liniar_space(next / 255.);
            double cur_val_gamma = to_liniar_space(cur_val / (num_of_colors - 1));
            double diff = next_gamma - prev_gamma;
            double dith_error = error_get_func(x, y);
            double new_value_with_error = cur_val_gamma + diff * dith_error;
            uchar result = (new_value_with_error >= next_gamma) ? next : prev;
            output.output << result;
        }
    }
}

void ordered_dithering() {

    static const std::vector<std::vector<double>> threshold_matrix = {
            {0,  32, 8,  40, 2,  34, 10, 42},
            {48, 16, 56, 24, 50, 18, 58, 26},
            {12, 44, 4,  36, 14, 46, 6,  38},
            {60, 28, 52, 20, 62, 30, 54, 22},
            {3,  35, 11, 43, 1,  33, 9,  41},
            {51, 19, 59, 27, 49, 17, 57, 25},
            {15, 47, 7,  39, 13, 45, 5,  37},
            {63, 31, 55, 23, 61, 29, 53, 21}
    };

    static const auto get_error = [](int x, int y) -> double {
        return threshold_matrix[y % 8][x % 8] / 64.;
    };
    dithering_error_not_saving(get_error);
}

void random_dithering() {
    std::mt19937 gen(time(0));
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    auto get_error = [&](int x, int y) {
        double d = dis(gen);
        while (d == 1.0) {
            d = dis(gen);
        }
        return d;
    };
    dithering_error_not_saving(get_error);
}

void dithering_error_saving(const std::vector<std::vector<int>> &parts, int divisor) {

    auto in_range = [](int x, int y) {
        if (x >= 0 && x < input.width) {
            if (y >= 0 && y < input.height) {
                return true;
            }
        }
        return false;
    };

    auto save_error = [&](int x_val, int y_val, double part, double quant_error) {
        if (in_range(x_val, y_val)) {
            error_buffer.get_ref(x_val, y_val) =
                    error_buffer.get_ref(x_val, y_val) + quant_error * part / divisor;
        }
    };

    for (int y = 0; y < input.height; y++) {
        for (int x = 0; x < input.width; ++x) {
            double cur_val = gradient_buffer.get_ref(x, y);
            auto[prev, next] = give_neighbours(cur_val);
            double prev_gamma = to_liniar_space(prev / 255.);
            double next_gamma = to_liniar_space(next / 255.);
            double cur_val_gamma = to_liniar_space(cur_val / (num_of_colors - 1));
            double cur_error = error_buffer.get_ref(x, y);
            double cur_value_with_error = cur_val_gamma + cur_error;
            double middle = (prev_gamma + next_gamma) / 2.;
            uchar result = prev;
            if (cur_value_with_error >= next_gamma) {
                result = next;
                cur_error = (cur_value_with_error - next_gamma);
            } else {
                cur_error = cur_value_with_error - prev_gamma;
            }
            output.output << result;
            for (int i = -2; i <= 2; i++) {
                for (int j = 0; j <= 2; ++j) {
                    save_error(x + i, y + j, parts[j][i + 2], cur_error);
                }
            }
        }
    }
}

void floyd_steinberg_dithering() {
    dithering_error_saving(
            {
                    {0, 0, 0, 7, 0},
                    {0, 3, 5, 1, 0},
                    {0, 0, 0, 0, 0},
            },
            16);
}

void jarvis_judice_ninke_dithering() {
    dithering_error_saving(
            {
                    {0, 0, 0, 7, 5},
                    {3, 5, 7, 5, 3},
                    {1, 3, 5, 3, 1},
            },
            48);

}

void sierra_dithering() {
    dithering_error_saving(
            {
                    {0, 0, 0, 5, 3},
                    {2, 4, 5, 4, 2},
                    {0, 2, 3, 2, 0},
            },
            32);
}

void atkinson_dithering() {
    dithering_error_saving(
            {
                    {0, 0, 0, 1, 1},
                    {0, 1, 1, 1, 0},
                    {0, 0, 1, 0, 0},
            },
            8);
}

void halftone_dithering() {
    static const std::vector<std::vector<double>> threshold_matrix = {
            {7,  13, 11, 4},
            {12, 16, 14, 8},
            {10, 15, 6,  2},
            {5,  9,  3,  1}
    };
    static const auto get_error = [](int x, int y) -> double {
        return threshold_matrix[y % 4][x % 4] / 16.;
    };
    dithering_error_not_saving(get_error);
}

using modife_func_ptr_t = void (*)();

const std::vector<modife_func_ptr_t> dithering_functions = {
        &no_dithering,
        &ordered_dithering,
        &random_dithering,
        &floyd_steinberg_dithering,
        &jarvis_judice_ninke_dithering,
        &sierra_dithering,
        &atkinson_dithering,
        &halftone_dithering
};

void fill_vector() {
    std::vector<uchar> nums(num_of_colors);
    num_to_byte.resize(num_of_colors, 0);
    for (int i = 0; i < num_of_colors; ++i) {
        nums[i] = (i << (8 - input_bits));
    }
    int start_pos = 7;
    for (; start_pos >= 0; start_pos -= input_bits) {
        for (int i = 0; i < num_of_colors; ++i) {
            num_to_byte[i] |= nums[i];
            nums[i] >>= input_bits;
        }
    }
}

bool read_arguments(int argc, char *argv[]) {
    static std::string usage = "Usage : lab3.exe <name_of_input_file> <name_of_output_files> <gradient> <dithering> <bits> <gamma>\n";
    if (argc != 7) {
        std::cerr << "Incorrect arguments number\n" << usage;
        return false;
    }
    input_file_name = argv[1];
    output_file_name = argv[2];
    try {
        input_gradient = std::stoi(argv[3]);
        if (input_gradient != 0 && input_gradient != 1) {
            std::cerr << "Incorrect number for gradient(value should be 0 or 1)\n" << usage;
            return false;
        }
        input_dithering = std::stoi(argv[4]);
        if (input_dithering < 0 || input_dithering > 7) {
            std::cerr << "Incorrect number for dithering(value should be in range [0..7])\n" << usage;
            return false;
        }
        input_bits = std::stoi(argv[5]);
        if (input_bits < 1 || input_bits > 8) {
            std::cerr << "Incorrect number for bits(value should be in range [1..8])\n" << usage;
            return false;
        }
        num_of_colors = (1 << input_bits);
        fill_vector();
        input_gamma = std::stod(argv[6]);
    } catch (std::invalid_argument &ex) {
        std::cerr << "Cannot parse argument as number\n" << usage;
        return false;
    }
    return true;
}

int main(int argc, char *argv[]) {
    bool testing = true;

    if (testing) {
        int argct = 7;
        char **argvt = new char *[argct];
        argvt[0] = "lab3.exe";
        argvt[1] = "B:\\Projects\\GitProjects\\Graphics\\pictures\\input_pictures\\field_6.pgm";
        argvt[2] = "B:\\Projects\\GitProjects\\Graphics\\pictures\\output_pictures\\field_2_6.pgm";
        argvt[3] = "1"; // gradient 0 or 1
        argvt[4] = "4"; // dithering
        argvt[5] = "3"; // bits 1..8
        argvt[6] = "2.2"; // 0 - sRGB, or else
        /*
0 – Нет дизеринга;
1 – Ordered (8x8);
2 – Random;
3 – Floyd–Steinberg;
4 – Jarvis, Judice, Ninke;
5 - Sierra (Sierra-3);
6 - Atkinson;
7 - Halftone (4x4, orthogonal);
         */
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
    auto &t = input;
    if (!input.read_header()) {
        return 1;
    }
    output.output_file_name = output_file_name;
    if (!output.open_output_file()) {
        return 1;
    }
    output.set_header_for_output_file(true, input.width, input.height, 255);
    output.write_header();
    if (!gradient_buffer.alloc(input.width, input.height)) {
        return 1;
    }
    if (input_dithering <= 6 && input_dithering >= 3) {
        if (!error_buffer.alloc(input.width, input.height)){
            return 1;
        }
        error_buffer.fill_zeros();
    }
    if (input_gradient == 1) {
        gradient_buffer.fill_gradient();
    } else {
        if (!input.read_buffer()){
            return 1;
        }
        gradient_buffer.fill_with_data();
    }
    dithering_functions[input_dithering]();
    return 0;
}
