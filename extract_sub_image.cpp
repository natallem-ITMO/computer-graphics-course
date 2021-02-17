#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

std::string format;
size_t width;
size_t height;
size_t max_color_size;

size_t x0;
size_t y0;
size_t sub_width;
size_t sub_height;

void check_size(std::ifstream &input);

bool read_arguments(int argc, char *argv[]);

void read_header(std::ifstream &buff);

std::vector<char> read_image(std::ifstream &input);

void write_header(std::ofstream &ofstream);

void write_sub_image(std::ofstream &output, std::vector<char> &vector);

size_t IX(size_t x, size_t y) {
    return x + y * width;
}

int main(int argc, char *argv[]) {
    if (!read_arguments(argc, argv)) {
        return 0;
    }

    std::ifstream input("..\\pictures\\hello.pgm");

    read_header(input);

    if (x0 + sub_width > width || y0 + sub_height > height) {
        std::cerr << "Sub image cannot be created because of size of sub image";
    } else {
        size_t size = width * height;
        std::vector<char> bytes = read_image(input);
        input.read(&bytes[0], size);
        std::ofstream output("..\\pictures\\hello_sub_image.pgm");
        write_header(output);
        write_sub_image(output, bytes);
        output.close();
    }
    input.close();
    return 0;
}

bool read_arguments(int argc, char *argv[]) {
    if (argc != 4) {
        std::cerr << "Incorrect arguments number\n";
        return false;
    }
    x0 = std::stoi(argv[1]);
    y0 = std::stoi(argv[2]);
    sub_width = std::stoi(argv[3]);
    sub_height = std::stoi(argv[4]);
    return true;
}

void read_header(std::ifstream &buff) {
    buff >> format;
    buff >> width;
    buff >> height;
    buff >> max_color_size;
    buff.get();
}

std::vector<char> read_image(std::ifstream &input) {
    size_t size = width * height;
    std::vector<char> bytes(size);
    input.read(&bytes[0], size);
    return bytes;
}

void write_header(std::ofstream &ofstream) {
    ofstream << "P5\n" + std::to_string(sub_width) + " " + std::to_string(sub_height) + "\n255\n";
}

void write_sub_image(std::ofstream &output, std::vector<char> &vector) {
    for (int j = y0; j < y0 + sub_height; ++j) {
        for (int i = x0; i < sub_width + x0; ++i) {
            output << (unsigned char) vector[IX(i, j)];
        }
    }
}

void check_size(std::ifstream &input) {
    size_t length = input.tellg();
    input.seekg(0, std::ios::end);
    size_t length1 = input.tellg();
    std::cout << length1 - length << "\n";
}

