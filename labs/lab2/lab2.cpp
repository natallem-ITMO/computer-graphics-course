#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <tuple>

int input_color_space = -1;
int output_color_space = -1;
int number_of_input_files = 0;
std::string input_file_name;
int number_of_output_files = 0;
std::string output_file_name;

bool is_P5 = true;
size_t width;
size_t height;
unsigned char max_color_size;
size_t size;

char *input_buffer_1 = nullptr;
char *input_buffer_2 = nullptr;
char *input_buffer_3 = nullptr;

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

double constraint(double d) {
//    if (d > 255)
//        return 255;
//    if (d < 0)
//        return 0;
    return d;
}

constexpr bool compressed = false;
struct RGB;
struct YCbCr_601;

struct RGB {
    double R = -1;
    double G = -1;
    double B = -1;

    void from_YCbCr_601(YCbCr_601 sp);

    [[nodiscard]] std::tuple<int, int, int> get_int_tuple() const;
};

std::ostream &operator<<(std::ostream &os, const RGB &dt) {
    os << "R=" << round(dt.R) << " G=" << round(dt.G) << " B=" << round(dt.B) << '\n';
    return os;
}

struct YCbCr_601 {
    double Y = -1;
    double Cb = -1;
    double Cr = -1;

    void from_RBG(RGB sp);

    YCbCr_601 compress();

    void decompress();

    std::tuple<int, int, int> get_int_tuple() const;

};

std::ostream &operator<<(std::ostream &os, const YCbCr_601 &dt) {
    os << "Y=" << round(dt.Y) << " Cb=" << round(dt.Cb) << " Cr=" << round(dt.Cr) << '\n';
    return os;
}

void RGB::from_YCbCr_601(YCbCr_601 input) {
    YCbCr_601 sp = input.compress();
    R = constraint(1.164 * (sp.Y - 16) + 1.596 * (sp.Cr - 128));
    G = constraint(1.164 * (sp.Y - 16) - 0.813 * (sp.Cr - 128) - 0.391 * (sp.Cb - 128));
    B = constraint(1.164 * (sp.Y - 16) + 2.018 * (sp.Cb - 128));
}

std::tuple<int, int, int> RGB::get_int_tuple() const {
    return std::make_tuple(round(R), round(G), round(B));
}

void YCbCr_601::from_RBG(RGB sp) {
    Y = constraint(0.257 * sp.R + 0.504 * sp.G + 0.098 * sp.B + 16);
    Cb = constraint(-0.148 * sp.R - 0.291 * sp.G + 0.439 * sp.B + 128);
    Cr = constraint(0.439 * sp.R - 0.368 * sp.G - 0.071 * sp.B + 128);
    if constexpr (compressed) {
        decompress();
    }
}

YCbCr_601 YCbCr_601::compress() {
    if constexpr (compressed) {
        return YCbCr_601{Y * (235 - 16) / 255. + 16, Cb * (240. - 16) / 255. + 16, Cr * (240. - 16) / 255. + 16};
    } else {
        return YCbCr_601(*this);
    }
}

void YCbCr_601::decompress() {
    Y = (Y - 16) / (235. - 16) * 255.;
    Cb = (Cb - 16) / (240. - 16) * 255.;
    Cr = (Cr - 16) / (240. - 16) * 255.;
}

std::tuple<int, int, int> YCbCr_601::get_int_tuple() const {
    return std::make_tuple(round(Y), round(Cb), round(Cr));
}

void check_RGB_init() {
    for (int r = 0; r <= 255; r++) {
        for (int g = 0; g <= 255; g++) {
            for (int b = 0; b <= 255; b++) {
                RGB red{r, g, b};
                YCbCr_601 d;
                d.from_RBG(red);
                red.from_YCbCr_601(d);
                if (red.get_int_tuple() != std::make_tuple(r, g, b)) {
                    std::cout << "not_equal:\n" << RGB{r, g, b} << red << "\n";
                }
            }
        }
    }
    RGB red{0,3,255};
    std::cout << "init " << red;
    YCbCr_601 d;
    d.from_RBG(red);
    std::cout << "transformed " << d;
    red.from_YCbCr_601(d);
    std::cout << "back result " << red << '\n';
}

void check_YCbCr_init() {
   /* for (int y = 0; y <= 255; y++) {
        for (int cb = 0; cb <= 255; cb++) {
            for (int cr = 0; cr <= 255; cr++) {
                YCbCr_601 init{y, cb, cr};
                RGB trans;
                trans.from_YCbCr_601(init);
                init.from_RBG(trans);
                if (init.get_int_tuple() != std::make_tuple(y, cb, cr)) {
                    std::cout << "not_equal:\n" << YCbCr_601{y, cb, cr} << init << "\n";
                }
            }
        }
    }*/
    YCbCr_601 d{235,240,220};
    std::cout << "init " << d;
    RGB red;
    red.from_YCbCr_601(d);
    std::cout << "transformed " << red;
    d.from_RBG(red);
    std::cout << "back result " << d << '\n';
}

int main(int argc, char *argv[]) {
//    check_RGB_init();
   check_YCbCr_init();
//601


    /* if (!read_arguments(argc, argv)) {
         return 1;
     }

     std::ifstream input;
     std::ofstream output;
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
         std::cerr << "Allocation for input image buffer failed: " << e.what() << '\n';
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
     return 0;*/
}

bool read_arguments(int argc, char *argv[]) {
    /*  static std::string usage = "Usage : lab2.exe -f <input_color_space> -t <output_color_space> "
                                 "-i <number_of_input_files> <name_of_input_file> "
                                 "-o <number_of_output_files> <name_of_output_files>\n";
      if (argc != 11) {
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
      return true;*/
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
    /*   ofstream << (is_P5 ? "P5" : "P6") << '\n';
       if (mode == 3 || mode == 4) {
           ofstream << std::to_string(height) + " " + std::to_string(width) << '\n';
       } else {
           ofstream << std::to_string(width) + " " + std::to_string(height) << '\n';
       }
       ofstream << std::to_string(max_color_size) << '\n';*/
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
                size_t cur_w = (width - 1 - w) * 3;
                ofstream << input_bytes[IX(cur_w, h)];
                ofstream << input_bytes[IX(cur_w + 1, h)];
                ofstream << input_bytes[IX(cur_w + 2, h)];
            }
        }
    }
}
