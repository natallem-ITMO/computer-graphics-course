#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <tuple>
#include <valarray>
#include <cassert>
#include <complex>
#include <string>

const bool need_to_transform = true; // for debug only to compare with OpenCV result
int Ncl;
int Mcl;

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

    unsigned char get_color(int x, int y) {
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

using base = std::complex<double>;
const double PI = 3.141592653589793238460;

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

void amplitude_phase_correlation(std::vector<std::vector<base>> &f_1, std::vector<std::vector<base>> &f_2,
                                 std::vector<std::vector<base>> &q_coefs) {
    for (int u = 0; u < Mcl; ++u) {
        for (int v = 0; v < Ncl; ++v) {
            q_coefs[u][v] = f_1[u][v] * std::conj(f_2[u][v]);
        }
        fft(q_coefs[u], false);
    }
}

void phase_correlation(std::vector<std::vector<base>> &f_1, std::vector<std::vector<base>> &f_2,
                       std::vector<std::vector<base>> &q_coefs) {
    for (int u = 0; u < Mcl; ++u) {
        for (int v = 0; v < Ncl; ++v) {
            q_coefs[u][v] = f_1[u][v] * std::conj(f_2[u][v]) / (std::abs(f_1[u][v]) * std::abs(std::conj(f_2[u][v])));
        }
        fft(q_coefs[u], false);
    }
}

using correlation_func_t = void (*)(std::vector<std::vector<base>> &, std::vector<std::vector<base>> &,
                                    std::vector<std::vector<base>> &);

std::vector<correlation_func_t> correlation_funcs = {&amplitude_phase_correlation, &phase_correlation};

void calc_correlation_fast() {
    int M = input_1.width;
    int N = input_1.height;

    std::vector<std::vector<base>> r_coefs_1(Ncl, std::vector<base>(Mcl, 0));
    std::vector<std::vector<base>> r_coefs_2(Ncl, std::vector<base>(Mcl, 0));
    for (int y = 0; y < Ncl; ++y) {
        if (y < N) {
            for (int x = 0; x < M; ++x) {
                r_coefs_1[y][x] = (double) input_1.get_color(x, y);
                r_coefs_2[y][x] = (double) input_2.get_color(x, y);
            }
        }
        fft(r_coefs_1[y], true);
        fft(r_coefs_2[y], true);
    }

    std::vector<std::vector<base>> f_values_1(Mcl, std::vector<base>(Ncl));
    std::vector<std::vector<base>> f_values_2(Mcl, std::vector<base>(Ncl));
    for (int u = 0; u < Mcl; u++) {
        for (int y = 0; y < Ncl; ++y) {
            f_values_1[u][y] = r_coefs_1[y][u];
            f_values_2[u][y] = r_coefs_2[y][u];
        }
        fft(f_values_1[u], true);
        fft(f_values_2[u], true);
    }

    std::vector<std::vector<base>> q_coefs(Mcl, std::vector<base>(Ncl));
    correlation_funcs[correlation_way](f_values_1, f_values_2, q_coefs);

    std::vector<std::vector<base>> Cf1f2(Ncl, std::vector<base>(Mcl, 0));
    for (int dy = 0; dy < Ncl; ++dy) {
        for (int u = 0; u < Mcl; ++u) {
            Cf1f2[dy][u] = q_coefs[u][dy];
        }
        fft(Cf1f2[dy], false);
    }
    auto mmax = Cf1f2[0][0].real();
    auto mmin = Cf1f2[0][0].real();
    std::complex<int> max_p(0, 0);

    for (int dy = 0; dy < Ncl; ++dy) {
        for (int dx = 0; dx < Mcl; ++dx) {
            double d = Cf1f2[dy][dx].real();
            if (d < mmin) {
                mmin = d;
            }
            if (d > mmax) {
                mmax = d;
                max_p = std::complex<int>(dx, dy);
            }
        }
    }

    int f = Mcl / 2;
    int g = Ncl / 2;
    std::cout << "(" << (max_p.real() + f) % Mcl << "," << (max_p.imag() + g) % Ncl << ")\n";
    for (int y = 0; y < Ncl; ++y) {
        for (int x = 0; x < Mcl; ++x) {
            if (need_to_transform) {
                double v = Cf1f2[(y + g) % Ncl][(x + f) % Mcl].real();
                double cur_part = (v - mmin) / (mmax - mmin);
                cur_part *= 255;
                output.output << (unsigned char) (cur_part);
            } else {
                double v = Cf1f2[y][x].real();
                double cur_part = (v - mmin) / (mmax - mmin);
                cur_part *= 255;
                output.output << (unsigned char) (cur_part);
            }
        }
    }

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
    if (!read_arguments(argc, argv)) {
        return 1;
    }
    input_1.input_file_name = input_file_name_1;
    if (!input_1.open_input_file() || !input_1.read_header() || !input_1.read_buffer()) {
        return 1;
    }
    input_2.input_file_name = input_file_name_2;
    if (!input_2.open_input_file() || !input_2.read_header() || !input_2.read_buffer()) {
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
    Mcl = find_closest(input_1.width);
    Ncl = find_closest(input_1.height);
    output.set_header_for_output_file(input_1.is_P5, Mcl, Ncl, 255);
    output.write_header();
    calc_correlation_fast();
    return 0;
}

