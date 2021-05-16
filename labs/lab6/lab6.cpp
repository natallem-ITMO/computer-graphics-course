#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <tuple>
#include <cassert>
#include <functional>
#include <complex>
#include <map>
#include <string>
#include <unordered_set>

static const double MY_PI = atan(1) * 4;

using uchar = unsigned char;

std::string input_file_name_1;
std::string input_file_name_2;
std::string output_file_name;
int correlation_way;

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

    int normalize(int x, int width) {
        if (x < 0) {
            x += width;
        }
        if (x >= width) {
            x -= width;
        }
        return x;
    }

    unsigned char get_color(int x, int y) {
        int xx = normalize(x, width);
        int yy = normalize(y, height);
        int wid = width;
        int he = height;
        if (!(xx >= 0 && xx < width && yy >= 0 && yy < height)) {
            int x = 10;
        }
        assert(xx >= 0 && xx < width && yy >= 0 && yy < height);
        return input_buffer[xx + yy * width];
    }

    unsigned char get_color_simple(int x, int y) {
        if (x >= width || y >= height) {
            x = x % width;
            y = y % height;
//            return 0;
        }
        return input_buffer[x + y * width];
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
file_worker input_1(program_error);
file_worker input_2(program_error);
file_worker output(program_error);

typedef std::complex<double> base;
static const double PI = atan(1) * 4;

void fft(std::vector<base> &a, bool invert) {
    int n = (int) a.size();
    if (n == 1) return;

    std::vector<base> a0(n / 2), a1(n / 2);
    for (int i = 0, j = 0; i < n; i += 2, ++j) {
        a0[j] = a[i];
        a1[j] = a[i + 1];
    }
    fft(a0, invert);
    fft(a1, invert);

    double ang = 2 * PI / n * (invert ? -1 : 1);
    base w(1), wn(cos(ang), sin(ang));
    for (int i = 0; i < n / 2; ++i) {
        a[i] = a0[i] + w * a1[i];
        a[i + n / 2] = a0[i] - w * a1[i];
        w *= wn;
    }
}

int find_closest(int x) {
    int cur = 1;
    while (cur < x) {
        cur <<= 1;
    }
    return cur;
}

void calc_correlation_fast() {
    int K = input_1.height;
    int M = input_1.width;
    int N = std::max(K, M);
    N = find_closest(N);

    std::vector<std::vector<base>> r_coefs_1;
    std::vector<std::vector<base>> r_coefs_2;
    for (int y = 0; y < N; ++y) {
        std::vector<base> cur_f_coeffs_1(N, 0);
        std::vector<base> cur_f_coeffs_2(N, 0);
        for (int x = 0; x < N; ++x) {
            cur_f_coeffs_1[x] = (double) input_1.get_color_simple(x, y);
            cur_f_coeffs_2[x] = (double) input_2.get_color_simple(x, y);
        }
        fft(cur_f_coeffs_1, true);
        fft(cur_f_coeffs_2, true);
        r_coefs_1.push_back(std::move(cur_f_coeffs_1));
        r_coefs_2.push_back(std::move(cur_f_coeffs_2));
    }
    std::vector<std::vector<base>> f_values_1;
    std::vector<std::vector<base>> f_values_2;
    for (int u = 0; u < N; u++) {
        std::vector<base> cur_r_coeff_1(N, 0);
        std::vector<base> cur_r_coeff_2(N, 0);
        for (int y = 0; y < N; ++y) {
            cur_r_coeff_1[y] = r_coefs_1[y][u];
            cur_r_coeff_2[y] = r_coefs_2[y][u];
        }
        fft(cur_r_coeff_1, true);
        fft(cur_r_coeff_2, true);
        f_values_1.push_back(std::move(cur_r_coeff_1));
        f_values_2.push_back(std::move(cur_r_coeff_2));
    }

    std::vector<std::vector<base>> q_coefs;
    for (int u = 0; u < N; ++u) {
        std::vector<base> cur_coeff_fs(N, 0);
        for (int v = 0; v < N; ++v) {
            cur_coeff_fs[v] = f_values_1[u][v] * std::conj(f_values_2[u][v]);
        }
        fft(cur_coeff_fs, false);
        q_coefs.push_back(std::move(cur_coeff_fs));
    }
    std::vector<std::vector<base>> Cf1f2;
    for (int dy = 0; dy < N; ++dy) {
        std::vector<base> cur_q_coefs(N, 0);
        for (int u = 0; u < N; ++u) {
            cur_q_coefs[u] = q_coefs[u][dy];
        }
        fft(cur_q_coefs, false);
        Cf1f2.push_back(std::move(cur_q_coefs));
    }
    auto mmax = Cf1f2[0][0].real();
    auto mmin = Cf1f2[0][0].real();
    std::complex<int> max_p(0, 0);
    std::complex<int> min_p(0, 0);

    for (int dy = 0; dy < K; ++dy) {
        for (int dx = 0; dx < M; ++dx) {
            double d = Cf1f2[dy][dx].real();
            if (d < mmin) {
                mmin = d;
                min_p = std::complex<int>(dx, dy);
            }
            if (d > mmax) {
                mmax = d;
                max_p = std::complex<int>(dx, dy);
            }
        }
    }

    for (int y = 0; y < K; ++y) {
        for (int x = 0; x < M; ++x) {
            int xx = x;
            int yy = y;
            xx += M / 2;
            yy+= K / 2;
            double v = Cf1f2[yy % K][xx % M].real();
            double cur_part = (v - mmin) / (mmax - mmin);
            cur_part *= 255;
            output.output << (unsigned char) (cur_part);

        }
    }
   /* for (int y = 0; y < input_1.height; y++) {
        for (int x = 0; x < input_1.width; ++x) {
            int xx = x;
            x += M / 2;
            int yy = y;
            yy += K / 2;
            double v = Cf1f2[yy  % K][xx % M].real();
            double cur_part = (v - mmin) / (mmax - mmin);
            cur_part *= 255;
            output.output << (unsigned char) (cur_part);
        }
    }*/
}


void calc_correlation() {
    int max_y_1 = 0;
    int max_y_2 = 0;
    int min_y_1 = 1000000000;
    int min_y_2 = 1000000000;
    int max_x_1 = 0;
    int max_x_2 = 0;
    int min_x_1 = 1000000000;
    int min_x_2 = 1000000000;
    for (int y = 0; y < input_1.height; y++) {
        for (int x = 0; x < input_1.width; ++x) {
            int color1 = input_1.get_color(x, y);
            if (color1 == 0) {
                max_x_1 = std::max(x, max_x_1);
                min_x_1 = std::min(x, min_x_1);
                max_y_1 = std::max(y, max_y_1);
                min_y_1 = std::min(y, min_y_1);
            }
            int color2 = input_2.get_color(x, y);
            if (color2 == 0) {
                max_x_2 = std::max(x, max_x_2);
                min_x_2 = std::min(x, min_x_2);
                max_y_2 = std::max(y, max_y_2);
                min_y_2 = std::min(y, min_y_2);
            }
        }
    }
    printf("max_x_1=%d min_x_1=%d max_y_1=%d min_y_1=%d max_x_2=%d min_x_2=%d max_y_2=%d min_y_2=%d\n", max_x_1,
           min_x_1, max_y_1, min_y_1, max_x_2, min_x_2, max_y_2, min_y_2);
    int center_x = input_1.width / 2;
    int center_y = input_1.height / 2;
    int xxxx = input_1.width;
    int yyyy = input_1.height;
    std::vector<std::vector<double>> values(input_1.height, std::vector<double>(input_1.width));
    std::unordered_set<double> diff_vals;
    file_worker output_1(program_error);
    output_1.output_file_name = "B:\\Projects\\GitProjects\\Graphics\\pictures\\output_pictures\\second_output.pgm";
    output_1.open_output_file();
    output_1.set_header_for_output_file(input_1.is_P5, input_1.width, input_1.height, 255);
    output_1.write_header();

    for (int y = 0; y < input_1.height; y++) {
        for (int x = 0; x < input_1.width; ++x) {
            if (x == 15 && y == 17) {
                int x = 10;
            }
            int delta_x = x - center_x;
            int delta_y = center_y - y;
            double cur_value = 0;
            bool correct = true;
            for (int j = 0; j < input_1.height; ++j) {
                for (int i = 0; i < input_1.width; ++i) {

                    unsigned char t = input_1.get_color(i, j);
                    unsigned char t1 = input_2.get_color(i - delta_x, j - delta_y);
                    /*   if (x == 15 && y == 17) {
                           output.output << t;
                           output_1.output << (unsigned char) t1;
                       }*/

                    if (t != t1) {
                        correct = false;
                    }
                    cur_value += (input_1.get_color(i, j)) * (input_2.get_color(i - delta_x, j - delta_y));
                }
            }

            if (correct) {
                int xwer = 10;
            }
//            cur_value /= input_1.width * input_1.height;
            values[y][x] = cur_value;
            if (x == y) {
                printf("x=%d y=%d val=%f\n", x, y, cur_value);
            }
        }
    }
    double mmax = values[0][0];
    double mmin = values[0][0];
    std::complex<int> max_p(0, 0);
    std::complex<int> min_p(0, 0);
    for (int y = 0; y < input_1.height; y++) {
        for (int x = 0; x < input_1.width; ++x) {
            double d = values[y][x];
            if (d < mmin) {
                mmin = d;
                min_p = std::complex<int>(x, y);
            }
            if (d > mmax) {
                mmax = d;
                max_p = std::complex<int>(x, y);
            }
        }
    }
    for (int y = 0; y < input_1.height; y++) {
        for (int x = 0; x < input_1.width; ++x) {
            double v = values[y][x];
            double cur_part = (v - mmin) / (mmax - mmin);
            cur_part *= 255;
            output.output << (unsigned char) (cur_part);
        }
    }
    int delta_x = max_p.real() - center_x;
    int delta_y = center_y - max_p.imag();
    std::cout << "SIZE " << input_1.width << " " << input_1.height << "\n";
    std::cout << "delta " << delta_x << " " << delta_y << "\n";
    std::cout << "RESULT " << max_p.real() << " " << max_p.imag() << "\n";
}

bool read_arguments(int argc, char *argv[]) {
    static std::string usage = "Usage : lab6.exe <name_of_input_file_1> <name_of_input_file_2> <name_of_output_files> <correlation_way>\n";
    if (argc != 5) {
        std::cerr << "Incorrect arguments number\n" << usage;
        return false;
    }
    input_file_name_1 = argv[1];
    input_file_name_2 = argv[2];
    output_file_name = argv[3];
    auto &d = input_file_name_1;
    auto &d1 = input_file_name_2;
    try {
        correlation_way = std::stoi(argv[4]);
        if (correlation_way != 0 && correlation_way != 1) {
            std::cerr << "Incorrect number for correlation way (value should be 0 or 1)\n" << usage;
            return false;
        }
    } catch (std::invalid_argument &ex) {
        std::cerr << "Cannot parse argument as number\n" << usage;
        return false;
    }
    return true;
}


int main(int argc, char *argv[]) {
    std::vector<std::tuple<std::string, std::string, std::string>> names = {
            {"chair_1.pgm", "chair_2.pgm", "chair.pgm"},
            {"circle1.pgm", "circle2.pgm", "circle.pgm"},
            {"chair_1.pgm", "chair_2.pgm", "chair_fft.pgm"},
            {"circle1.pgm", "circle2.pgm", "circle_fft.pgm"},
//            {"circle1.pgm", "circle2.pgm", "circle_fft_origin.pgm"},
    };
    bool testing = true;
//    int i = names.size() - 1;
    int i = 2;
    double coef_w = 3;
    double coef_h = 4;
    if (testing) {
        std::string name_1 = std::get<0>(names[i]);
        std::string name_2 = std::get<1>(names[i]);
        std::string result = std::get<2>(names[i]);
        int argct = 5;
        char **argvt = new char *[argct];
        argvt[0] = "lab3.exe";
        std::string name_file = "B:\\Projects\\GitProjects\\Graphics\\pictures\\input_pictures\\" + name_1;
        std::string name_file_2 = "B:\\Projects\\GitProjects\\Graphics\\pictures\\input_pictures\\" + name_2;
        std::string name_file_res = "B:\\Projects\\GitProjects\\Graphics\\pictures\\output_pictures\\" + result;
        argvt[1] = const_cast<char *>(name_file.c_str());
        argvt[2] = const_cast<char *>(name_file_2.c_str());
        argvt[3] = const_cast<char *>(name_file_res.c_str());
        argvt[4] = "0";
        if (!read_arguments(argct, argvt)) {
            return 1;
        }
    } else {
        if (!read_arguments(argc, argv)) {
            return 1;
        }
    }
    input_1.input_file_name = input_file_name_1;
    if (!input_1.open_input_file()) {
        return 1;
    }
    if (!input_1.read_header()) {
        return 1;
    }
    if (!input_1.read_buffer()) {
        return 1;
    }
    input_2.input_file_name = input_file_name_2;
    if (!input_2.open_input_file()) {
        return 1;
    }
    if (!input_2.read_header()) {
        return 1;
    }
    if (!input_2.read_buffer()) {
        return 1;
    }
    if (!input_1.compare_with_other_input_file(input_2)) {
        std::cerr << "Two images for compare aren't equal in size or format\n";
        return 1;
    }
    output.output_file_name = output_file_name;
    if (!output.open_output_file()) {
        return 1;
    }
    output.set_header_for_output_file(input_1.is_P5, input_1.width, input_1.height, 255);
    output.write_header();
    calc_correlation_fast();
//    calc_correlation();
    return 0;
}