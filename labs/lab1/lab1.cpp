#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

std::string input_file_name;
std::string output_file_name;
int mode = 0;


bool is_P5 = true;
size_t width;
size_t height;
unsigned char max_color_size;
size_t size;

using modife_func_ptr_t = void (*)(std::ofstream &ofstream, char *input_bytes);

void inversion_modification(std::ofstream &ofstream, char *input_bytes);

void horizontal_flip_modification(std::ofstream &ofstream, char *input_bytes);

void vertical_flip_modification(std::ofstream &ofstream, char *input_bytes);

void rotate_counterclockwise_modification(std::ofstream &ofstream, char *input_bytes);

void rotate_clockwise_modification(std::ofstream &ofstream, char *input_bytes);

const std::vector<modife_func_ptr_t> modification_functions = {
        &inversion_modification,
        &horizontal_flip_modification,
        &vertical_flip_modification,
        &rotate_clockwise_modification,
        &rotate_counterclockwise_modification
};

bool read_arguments(int argc, char *argv[]);

bool read_header(std::ifstream &buff);

void write_header(std::ofstream &ofstream);

inline size_t IX(size_t x, size_t y) {
    return x + y * width * (is_P5 ? 1 : 3);
}

int main(int argc, char *argv[]) {
    if (!read_arguments(argc, argv)) {
        return 1;
    }

    std::ifstream input;
    std::ofstream output;
    char *input_image_buffer = nullptr;
    bool error = false;

    input.open(input_file_name, std::ios::binary);
    if (!input.is_open()) {
        std::cerr << "Failed to open input file \"" << input_file_name << "\"\n";
        error = true;
        goto free_resources;
    }

    if (!read_header(input)) {
        error = true;
        goto free_resources;
    }

    try {
        size = width * height * (is_P5 ? 1 : 3);
        input_image_buffer = new char[size];
        input.read(input_image_buffer, size);
        if (input.gcount() != size) {
            std::cerr << "Couldn't read all bytes of input image to buffer\n";
            error = true;
            goto free_resources;
        }
        output.open(output_file_name, std::ios::binary);
        if (!output.is_open()) {
            std::cerr << "Failed to open output file \"" << output_file_name << "\"\n";
            error = true;
            goto free_resources;
        }
        write_header(output);
        modification_functions[mode](output, input_image_buffer);
    } catch (const std::bad_alloc &e) {
        std::cerr << "Allocation for input image bugger failed: " << e.what() << '\n';
        error = true;
        goto free_resources;
    }

    free_resources:
    delete[] input_image_buffer;
    if (input.is_open()) {
        input.close();
    }
    if (output.is_open()) {
        output.close();
    }
    if (error) {
        return 1;
    }
    return 0;
}

bool read_arguments(int argc, char *argv[]) {
    static std::string usage = "Usage : lab1.exe <input_file_name> <output_file_name> <mode(integer number 0..4)>\n";
    if (argc != 4) {
        std::cerr << "Incorrect arguments number\n" << usage;
        return false;
    }
    try {
        input_file_name = argv[1];
        output_file_name = argv[2];
        mode = std::stoi(argv[3]);
    } catch (std::invalid_argument &ex) {
        std::cerr << "Cannot parse argument as integer\n" << usage;
        return false;
    }
    if (mode > 4 || mode < 0) {
        std::cerr << "Incorrect modification mode " << mode << '\n'
                  << "Mode can be only one of {0, 1, 2, 3, 4} numbers\n"
                  << usage;
        return false;
    }
    return true;
}

bool read_header(std::ifstream &buff) {
    std::string format;
    buff >> format;
    if (format != "P5" && format != "P6") {
        std::cerr << "Input image in file \"" << input_file_name << "\" is not in format P5 or P6\n";
        return false;
    }
    if (format != "P5") {
        is_P5 = false;
    }

    buff >> width;
    if (width == 0) {
        std::cerr << "Width of input image in file \"" << input_file_name << "\" is incorrect (= 0)\n";
        return false;
    }
    buff >> height;
    if (height == 0) {
        std::cerr << "Height of input image in file \"" << input_file_name << "\" is incorrect (= 0)\n";
        return false;
    }
    int temp;
    buff >> temp;
    max_color_size = temp;
    if (temp <= 0) {
        std::cerr << "Maximum color number of input image in file \"" << input_file_name << "\" is incorrect (<= 0)\n";
        return false;
    }
    buff.get();
    return true;
}

void write_header(std::ofstream &ofstream) {
    ofstream << (is_P5 ? "P5" : "P6") << '\n';
    if (mode == 3 || mode == 4) {
        ofstream << std::to_string(height) + " " + std::to_string(width) << '\n';
    } else {
        ofstream << std::to_string(width) + " " + std::to_string(height) << '\n';
    }
    ofstream << std::to_string(max_color_size) << '\n';
}

void inversion_modification(std::ofstream &ofstream, char *input_bytes) {
    for (size_t i = 0; i < size; i++) {
        ofstream << (unsigned char) (max_color_size - (unsigned char) input_bytes[i]);
    }
}

void horizontal_flip_modification(std::ofstream &ofstream, char *input_bytes) {
    if (is_P5) {
        for (size_t h = 0; h < height; h++) {
            for (size_t w = 0; w < width; w++) {
                ofstream << input_bytes[IX((width - 1 - w), h)];
            }
        }
    } else {
        for (size_t h = 0; h < height; h++) {
            for (size_t w = 0; w < width; w++) {
                size_t cur_w = (width - 1 - w) * 3;
                ofstream << input_bytes[IX(cur_w, h)];
                ofstream << input_bytes[IX(cur_w + 1, h)];
                ofstream << input_bytes[IX(cur_w + 2, h)];
            }
        }
    }
}

void vertical_flip_modification(std::ofstream &ofstream, char *input_bytes) {
    if (is_P5) {
        for (size_t h = 0; h < height; h++) {
            for (size_t w = 0; w < width; w++) {
                ofstream << input_bytes[IX(w, (height - 1 - h))];
            }
        }
    } else {
        for (size_t h = 0; h < height; h++) {
            for (size_t w = 0; w < width; w++) {
                size_t cur_h = height - 1 - h;
                size_t cur_w = w * 3;
                ofstream << input_bytes[IX(cur_w, cur_h)];
                ofstream << input_bytes[IX(cur_w + 1, cur_h)];
                ofstream << input_bytes[IX(cur_w + 2, cur_h)];
            }
        }
    }
}

void rotate_clockwise_modification(std::ofstream &ofstream, char *input_bytes) {
    if (is_P5) {
        for (size_t w = 0; w < width; w++) {
            for (size_t h = 0; h < height; h++) {
                ofstream << input_bytes[IX(w, (height - 1 - h))];
            }
        }
    } else {
        for (size_t w = 0; w < width; w++) {
            for (size_t h = 0; h < height; h++) {
                size_t cur_w = w * 3;
                size_t cur_h = height - 1 - h;
                ofstream << input_bytes[IX(cur_w, cur_h)];
                ofstream << input_bytes[IX(cur_w + 1, cur_h)];
                ofstream << input_bytes[IX(cur_w + 2, cur_h)];
            }
        }
    }
}

void rotate_counterclockwise_modification(std::ofstream &ofstream, char *input_bytes) {
    if (is_P5) {
        for (size_t w = 0; w < width; w++) {
            for (size_t h = 0; h < height; h++) {
                ofstream << input_bytes[IX(width - 1 - w, h)];
            }
        }
    } else {
        for (size_t w = 0; w < width; w++) {
            for (size_t h = 0; h < height; h++) {
                size_t cur_w = (width-1-w) * 3;
                ofstream << input_bytes[IX(cur_w, h)];
                ofstream << input_bytes[IX(cur_w + 1, h)];
                ofstream << input_bytes[IX(cur_w + 2, h)];
            }
        }
    }
}
