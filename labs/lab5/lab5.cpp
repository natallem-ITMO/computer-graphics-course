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

static const double MY_PI = atan(1) * 4;

using uchar = unsigned char;

std::string input_file_name;
std::string output_file_name;
int threshold_number;

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


bool read_arguments(int argc, char *argv[]) {
    static std::string usage = "Usage : lab3.exe <name_of_input_file> <name_of_output_files> <threshold_number>\n";
    if (argc != 4) {
        std::cerr << "Incorrect arguments number\n" << usage;
        return false;
    }
    input_file_name = argv[1];
    output_file_name = argv[2];
    try {
        threshold_number = std::stoi(argv[3]);
        if (threshold_number < 2 || threshold_number > 255) {
            std::cerr << "Incorrect number for threshold (value should be in range [1..255])\n" << usage;
            return false;
        }
        --threshold_number;
    } catch (std::invalid_argument &ex) {
        std::cerr << "Cannot parse argument as number\n" << usage;
        return false;
    }
    return true;
}

const static size_t color_number = 256;
std::vector<size_t> histogram(color_number, 0);
std::vector<double> possibilities(color_number);
std::vector<double> possibilities_sum(color_number);
std::vector<double> partial_expectation(color_number);
std::vector<double> partial_expectation_sum(color_number);
std::vector<int> thresholds;
std::vector<int> best_thresholds;
std::vector<int> best_thresholds_stupid;
double best_dispersion = 0;
double best_dispersion_stupid = 0;

void create_histogram() {
    for (int i = 0; i < input.number_of_pixels; ++i) {
        unsigned char val = input.input_buffer[i];
        ++histogram[val];
    }
}

void precalc() {
    double N = input.number_of_pixels;
    for (int i = 0; i < histogram.size(); ++i) {
        possibilities[i] = histogram[i] / N;
    }
    possibilities_sum[0] = possibilities[0];
    for (int i = 1; i < histogram.size(); ++i) {
        possibilities_sum[i] = possibilities_sum[i - 1] + possibilities[i];
    }
    partial_expectation[0] = 0;
    for (int i = 0; i < histogram.size(); ++i) {
        partial_expectation[i] = i * possibilities[i];
    }
    partial_expectation_sum[0] = partial_expectation[0];
    for (int i = 1; i < histogram.size(); ++i) {
        partial_expectation_sum[i] = partial_expectation_sum[i - 1] + partial_expectation[i];
    }
}

static const int L = 255;

double get_q_k(int k) {
    double first_sum = possibilities_sum[thresholds[k]];
    if (k == 1) {
        return first_sum;
    }
    double second_sum = possibilities_sum[thresholds[k - 1]];
//    return first_sum - second_sum;
    int start = thresholds[k - 1] + 1;
    int end = thresholds[k];
    double cur_sum = 0;
    for (int i = start; i <= end; ++i) {
        cur_sum += possibilities[i];
    }
    return first_sum - second_sum;
}

double get_q_k_stupid(int k) {
    int start = thresholds[k - 1] + 1;
    int end = thresholds[k];
    double cur_sum = 0;
    for (int i = start; i <= end; ++i) {
        cur_sum += possibilities[i];
    }
    return cur_sum;
}

double get_mu_k(int k) {
    double first_sum = partial_expectation_sum[thresholds[k]];
    if (k == 1) {
        return first_sum;
    }
    double second_sum = partial_expectation_sum[thresholds[k - 1]];
    return first_sum - second_sum;
}

double get_mu_k_stupid(int k) {
    int start = thresholds[k - 1] + 1;
    int end = thresholds[k];
    double cur_sum = 0;
    for (int i = start; i <= end; ++i) {
        cur_sum += partial_expectation[i];
    }
    return cur_sum;
}

double get_k_component(int k) {
    double q_k = get_q_k(k);
    double mu_k = get_mu_k(k);
    if (q_k == 0) {
        return 0;
    }
    return mu_k * mu_k / q_k;
}


double get_k_component_stupid(int k) {
    double q_k = get_q_k_stupid(k);
    double mu_k = get_mu_k_stupid(k);
    if (q_k == 0) {
        return 0;
    }
    return mu_k * mu_k / q_k;
}

template<typename T>
void show_vec(std::vector<T> &v, const std::string &message) {
    std::cout << message << ": ";
    for (auto t : v) {
        std::cout << t << ' ';
    }
    std::cout << "\n";
}


void check_cur_thresholds_stupid() {
//    show_vec(thresholds, "current thresholds");
    double cur_res = 0;
    for (int i = 1; i <= threshold_number + 1; ++i) {
        double d = get_k_component_stupid(i);
        cur_res += d;
    }
    if (cur_res > best_dispersion_stupid) {
        best_dispersion_stupid = cur_res;
        best_thresholds_stupid = thresholds;
    }
}


void check_cur_thresholds() {
//    show_vec(thresholds, "current thresholds");
    double cur_res = 0;
    for (int i = 1; i <= threshold_number + 1; ++i) {
        double d = get_k_component(i);
        cur_res += d;
/*        double d2 = get_k_component_stupid(i);
        double eps = 1e-5;
        std::cout << i << " component=" << d << " res=" << cur_res << " component_stupid=" << d2 << "\n";
        if (abs(d - d2) > eps){
            double  diff = abs(d-d2);
            int x =10;
            double d = get_k_component(i);
            double d2 = get_k_component_stupid(i);
        }*/
    }
    if (cur_res > best_dispersion) {
        best_dispersion = cur_res;
        best_thresholds = thresholds;
    }
}


void rec_threshold_calc(int cur_threshold_num, int prev_value) {
    for (int i = prev_value + 1; i <= L - threshold_number + (cur_threshold_num - 1); i++) {
        thresholds[cur_threshold_num] = i;
        if (cur_threshold_num != threshold_number) {
            rec_threshold_calc(cur_threshold_num + 1, i);
        } else {
            check_cur_thresholds();
            check_cur_thresholds_stupid();
        }
    }
    std::vector<int> &_best_thresholds = best_thresholds;
    std::vector<int> &_best_thresholds_stupid = best_thresholds_stupid;
    assert(best_thresholds == best_thresholds_stupid);
}

void find_best_thresholds() {
    std::vector<size_t> &_histogram = histogram;
    std::vector<double> &_possibilities = possibilities;
    std::vector<double> &_partial_expectation = partial_expectation;
    std::vector<int> &_thresholds = thresholds;
    std::vector<int> &_best_thresholds = best_thresholds;
    std::vector<int> &_best_thresholds_stupid = best_thresholds_stupid;
    double &_best_dispersion = best_dispersion;
    double &_best_dispersion_stupid = best_dispersion_stupid;
    create_histogram();
    precalc();
    thresholds.resize(threshold_number + 2);
    thresholds[0] = -1;
    thresholds[threshold_number + 1] = 255;
    rec_threshold_calc(1, -1);
}

void write_result() {
    std::map<unsigned char, unsigned char> color_map;
    color_map[best_thresholds[1]] = 0;
    color_map[best_thresholds[threshold_number + 1]] = 255;
    auto to_add = (unsigned char) (255. / threshold_number);
    unsigned char cur_value = to_add;
    for (int i = 2; i < threshold_number + 1; ++i) {
        color_map[best_thresholds[i]] = cur_value;
        cur_value += to_add;
    }
    for (size_t i = 0; i < input.number_of_pixels; ++i) {
        unsigned char cur_value = input.input_buffer[i];
        auto itr = color_map.lower_bound(cur_value);
        output.output << itr->second;
    }
}

int main(int argc, char *argv[]) {


    std::vector<std::pair<std::string, int>> names = {
            {"seeds.pgm",       2},
            {"big_seeds.pgm",   3},
            {"big_seeds.pgm",   2},
            {"bubbles.pgm",     3},
            {"bubbles.pgm",     2},
            {"black_faces.pgm", 2},
            {"black_faces.pgm", 3},
            {"woman_hat.pgm",   3},
            {"moon.pgm",        3},
            {"portal_tech.pgm", 3},
            {"road_woman.pgm",  3},
            {"road_woman.pgm",  4},
            {"image_1.pgm",     2},
            {"image_2.pgm",     2},
            {"image_3.pgm",     2},
            {"image_4.pgm",     2},
            {"woman_hat.pgm",   2},
//            {"woman_hat.pgm",   4},
    };
    bool testing = true;
    int i = names.size() - 1;
    if (testing) {
        bool examples = false;
        bool noised = false;
        std::string path_prefix;
        if (examples) {
            path_prefix += "7z_examples\\";
            if (noised) {
                path_prefix += "noised\\";
            } else {
                path_prefix += "blurred\\";
            }
        }


        std::string name = names[i].first;
        int classes = names[i].second;
        int argct = 4;
        char **argvt = new char *[argct];
        argvt[0] = "lab5.exe";
        std::string name_file = "B:\\Projects\\GitProjects\\Graphics\\pictures\\input_pictures\\" + path_prefix + name;
        argvt[1] = const_cast<char *>(name_file.c_str());
        argvt[3] = const_cast<char *>(std::to_string(classes).c_str()); // thresholds number
        std::string name_file_out =
                "B:\\Projects\\GitProjects\\Graphics\\pictures\\output_pictures\\lab_5_ex\\" + path_prefix +
                std::string(argvt[3]) + "_" +
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
    output.set_header_for_output_file(input.is_P5, input.width, input.height, 255);
    output.write_header();
    find_best_thresholds();
    write_result();
    return 0;
}
