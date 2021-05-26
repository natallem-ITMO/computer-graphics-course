#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <tuple>
#include <cassert>
#include <complex>
#include <string>
#include <map>
#include <set>
#include <memory>
#include <functional>
#include "zlib.h"

static int max_w = 40;
static int max_h = 20;
static int cur_w = 0;
static int cur_h = 0;

std::string input_file_name;
std::string output_file_name;
using uchar = unsigned char;

int byte_in_pixel;

struct file_worker {
    explicit file_worker() = default;

    bool open_input_file() {
        assert(!input_file_name.empty());
        input.open(input_file_name, std::ios::binary);
        if (!input.is_open()) {
            std::cerr << "Failed to open input file \"" << input_file_name << "\"\n";
            return false;
        }
        return true;
    }

    bool open_output_file() {
        output.open(output_file_name, std::ios::binary);
        if (!output.is_open()) {
            std::cerr << "Failed to open output file \"" << output_file_name << "\"\n";
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
            return false;
        }
        if (format != "P5") {
            return false;
        }
        input >> width;
        if (width == 0) {
            std::cerr << "Width of input image in file \"" << input_file_name << "\" is incorrect (= 0)\n";
            return false;
        }
        input >> height;
        if (height == 0) {
            std::cerr << "Height of input image in file \"" << input_file_name << "\" is incorrect (= 0)\n";
            return false;
        }
        int temp;
        input >> temp;
        max_color_size = temp;
        if (temp <= 0) {
            std::cerr << "Maximum color number of input image in file \"" << input_file_name
                      << "\" is incorrect (<= 0)\n";
            return false;
        }
        if (max_color_size != 255) {
            std::cerr << "Maximum color number of input image in file \"" << input_file_name
                      << "\" should be 255\n";
            return false;
        }
        input.get();
        size = width * height;
        number_of_pixels = (is_P5) ? size : size * 3;
        return true;
    }

    bool read_buffer() {
        try {
            input_buffer = new char[number_of_pixels];
            input.read(input_buffer, number_of_pixels);
            if (input.gcount() != number_of_pixels) {
                std::cerr << "Couldn't read all bytes of input image from file \'" << input_file_name
                          << "\' to buffer\n";
                return false;
            }
        } catch (const std::bad_alloc &e) {
            std::cerr << "Allocation for input image from file \'" << input_file_name << "\' buffer failed: "
                      << e.what() << '\n';
            return false;
        }
        return true;
    }

    uchar get_color(int x, int y) {
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

file_worker output;

struct chunk {
    chunk() = default;

    void swap(chunk &other) {
        std::swap(other.length, length);
        std::swap(other.type, type);
        std::swap(other.data, data);
        std::swap(other.CRC, CRC);
        std::swap(other.pos, pos);
    }

    chunk(chunk &&other) {
        length = other.length;
        type = std::move(other.type);
        data = other.data;
        other.data = nullptr;
        for (int i = 0; i < 4; ++i) {
            CRC[i] = other.CRC[i];
        }
        pos = other.pos;
    }

    bool read(std::ifstream &fin) {
        for (int i = 3; i >= 0; --i) {
            length = length | (fin.get() << (8 * i));
        }
        for (int i = 0; i < 4; i++) {
            char ch = fin.get();
            type[i] = ch;
        }
        try {
            data = new char[length];
            fin.read(data, length);
            int c = fin.gcount();
            if (c != length) {
                std::cerr << "Couldn't read all bytes of chunk from file \'" << input_file_name << "\' to buffer\n";
                return false;
            }
        } catch (const std::bad_alloc &e) {
            std::cerr << "Allocation for chunk data buffer failed: " << e.what() << '\n';
            return false;
        }
        for (int i = 0; i < 4; i++) {
            CRC[i] = fin.get();
        }
        return true;
    }

    ~chunk() {
        delete[] data;
    }

    bool get_bit() {
        if (pos_in_byte == 8) {
            cur_byte = get_byte();
            pos_in_byte = 0;
        }
        return cur_byte & (1 << (pos_in_byte++));
    }

    char get_byte() {
        assert(pos_in_byte == 8);
        if (pos >= length) {
            throw std::out_of_range("no readable elements in chunk");
        }
        return data[pos++];
    }

    int get_n_bits_reversed_order(int n) {
        int res = 0;
        std::string s;
        for (int i = 0; i < n; i++) {
            bool ch = get_bit();
            s += (ch) ? '1' : '0';
            res |= (ch << i);
        }
        return res;
    }

    int get_n_bits_direct_order(int n) {
        int res = 0;
        std::string s;
        for (int i = 0; i < n; i++) {
            bool chh = get_bit();
            res <<= 1;
            res |= chh;
            s += (chh) ? '1' : '0';
        }
        return res;
    }


    uint32_t length = 0;
    std::string type = "1234";
    char *data = nullptr;
    uchar CRC[4];
    size_t pos = 0;
    char cur_byte;
    int pos_in_byte = 8;
};

struct zlib_decoder {

    zlib_decoder(chunk &ch, char *&outputBuffer, size_t &outputBufferPos, size_t max_buffer_size)
            : ch(ch), output_buffer(outputBuffer), output_buffer_pos(outputBufferPos),
              max_buffer_size(max_buffer_size) {}

    bool read_zlib_header() {
        int CMF = ch.get_byte();
        CM = 0b00001111 & CMF;
        if (CM != 8) {
            std::cerr << "CM in zlib header != 8\n";
            return false;
        }
        unsigned char FLG = ch.get_byte();
        FLEVEL = FLG & 0b11000000;
        FDICT = FLG & 0b00100000;
        FCHECK = FLG & 0b00011111;
        int diff = FLG - FCHECK;
        int not_divable_31 = ((256 * CMF) + diff);
        assert((not_divable_31 % 31) != 0);
        int x = not_divable_31 / 31;
        ++x;
        int new_diff = (x * 31) - not_divable_31;
        assert(new_diff == FCHECK);
        return true;
    }

    using decode_func_t = bool (zlib_decoder::*)();

    bool decode_data() {
        bool end_of_blocks = false;
        while (!end_of_blocks) {
            end_of_blocks = ch.get_bit();
            int mode = ch.get_n_bits_reversed_order(2);
            if (mode >= 3) {
                std::cerr << "Unsupported mode for zlib block = " << mode << "\n";
                return false;
            }
            if (!(this->*decode_funcions[mode])()) {
                return false;
            }
        }
        return true;
    }

    bool decode_no_compression() {
        ch.pos_in_byte = 8;
        char b1 = ch.get_byte();
        char b2 = ch.get_byte();
        int length = (b1 << 8) | b2;
        char b1c = ch.get_byte();
        char b2c = ch.get_byte();
        assert((b1c & b1) == 0b0000'0000'0000'0000);
        assert((b2c & b2) == 0b0000'0000'0000'0000);
        assert((b1c | b1) == 0b1111'1111'1111'1111);
        assert((b2c | b2) == 0b1111'1111'1111'1111);
        for (int i = 0; i < length; ++i) {
            write_to_buffer(ch.get_byte());
        }
        return true;
    }

    bool decode_Huffman_codes(const std::function<int()> &get_command, const std::function<int()> &get_shift) {
        static const std::vector<std::pair<int, int>> length_code_table = {
                {0, 3},
                {0, 4},
                {0, 5},
                {0, 6},
                {0, 7},
                {0, 8},
                {0, 9},
                {0, 10},
                {1, 11},
                {1, 13},
                {1, 15},
                {1, 17},
                {2, 19},
                {2, 23},
                {2, 27},
                {2, 31},
                {3, 35},
                {3, 43},
                {3, 51},
                {3, 59},
                {4, 67},
                {4, 83},
                {4, 99},
                {4, 115},
                {5, 131},
                {5, 163},
                {5, 195},
                {5, 227},
                {0, 258}
        };
        static const std::vector<std::pair<int, int>> dist_code_table = {
                {0,  1},
                {0,  2},
                {0,  3},
                {0,  4},
                {1,  5},
                {1,  7},
                {2,  9},
                {2,  13},
                {3,  17},
                {3,  25},
                {4,  33},
                {4,  49},
                {5,  65},
                {5,  97},
                {6,  129},
                {6,  193},
                {7,  257},
                {7,  385},
                {8,  513},
                {8,  769},
                {9,  1025},
                {9,  1537},
                {10, 2049},
                {10, 3073},
                {11, 4097},
                {11, 6145},
                {12, 8193},
                {12, 12289},
                {13, 16385},
                {13, 24577},
        };
        std::string debug_res;
        while (true) {
            int actual_value = get_command();
            if (actual_value == 256) {
                break;
            }
            assert(actual_value != 286 || actual_value != 287);
            if (actual_value > 256) {
                int length = length_code_table[actual_value - 257].second;
                int add_bits = length_code_table[actual_value - 257].first;
                if (add_bits != 0) {
                    length += ch.get_n_bits_reversed_order(add_bits);
                }
                int table_shift = get_shift();
                int shift = dist_code_table[table_shift].second;
                add_bits = dist_code_table[table_shift].first;
                int prev_shift = shift;
                int k;
                if (add_bits != 0) {
                    k = ch.get_n_bits_reversed_order(add_bits);
                    shift += k;
                }
//                std::cout << " length = " << length << " dist = " << shift;
                size_t initial_buffer_pos = output_buffer_pos;
                size_t shifted_position = initial_buffer_pos - shift;
                for (int j = 0; j < length; j++) {
                    if (shifted_position == initial_buffer_pos) {
                        shifted_position = initial_buffer_pos - shift;
                    }
                    if (cur_w == 40 && cur_h == 19) {
                        int x = 10;
                    }
                    if (output_buffer_pos >= max_buffer_size) {
                        break;
                    }
                    write_to_buffer(output_buffer[shifted_position++]);
                }
            } else {
                write_to_buffer((uchar) actual_value);
            }
        }
        return true;
    }

    bool decode_fixed_Huffman_codes() {

        const static std::vector<std::tuple<int, int, std::pair<int, int>>> code_table = {
                {256, 7, {0b0000000,   0b0010111}},
                {0,   8, {0b00110000,  0b10111111}},
                {280, 8, {0b11000000,  0b11000111}},
                {144, 9, {0b110010000, 0b111111111}},
        };

        auto get_command_func = [this]() -> int {
            int cur = ch.get_n_bits_direct_order(7);
            int i;
            for (i = 0; i < 4; ++i) {
                int k = std::get<2>(code_table[i]).first;
                int d = std::get<2>(code_table[i]).second;
                if (cur >= std::get<2>(code_table[i]).first && cur <= std::get<2>(code_table[i]).second) {
                    break;
                } else if (i != 1) {
                    cur <<= 1;
                    cur |= ch.get_bit();
                }
            }
            return cur - std::get<2>(code_table[i]).first + std::get<0>(code_table[i]);
        };

        auto get_dist_func = [this]() -> int {
            return ch.get_n_bits_direct_order(5);
        };
        return decode_Huffman_codes(get_command_func, get_dist_func);
    }

    struct node {
        int val = -1;
        node *left = nullptr;
        node *right = nullptr;

        ~node() {
            delete left;
            delete right;
        }
    };

    void create_branch(node *cur_node, int value, int pos, std::vector<bool> &cur_num) {
        if (pos == cur_num.size()) {
            assert(cur_node->val == -1);
            cur_node->val = value;
            return;
        }
        if (cur_num[pos]) {
            if (!cur_node->right) {
                cur_node->right = new node();
            }
            create_branch(cur_node->right, value, pos + 1, cur_num);
        } else {
            if (!cur_node->left) {
                cur_node->left = new node();
            }
            create_branch(cur_node->left, value, pos + 1, cur_num);
        }
    }

    std::unique_ptr<node> create_Huffman_tree(std::map<int, std::set<int>> &mapa) {
        mapa.erase(0);
        if (mapa.size() == 1 && (*mapa.begin()).second.size() == 1) {
            int x = 10;
        }
        assert(mapa.size() != 1 || (*mapa.begin()).second.size() != 1);
        std::map<int, std::vector<bool>> values;
        int cur_value = 0;
        int prev_len = 0;
        for (auto &p : mapa) {
            int len = p.first;
            int diff_len = len - prev_len;
            cur_value <<= diff_len;
            for (auto val : p.second) {
                values[val].resize(len);
                for (int i = 0; i < len; ++i) {
                    values[val][i] = (cur_value & (1 << (len - 1 - i))) != 0;
                }
                ++cur_value;
            }
            prev_len = len;
        }
        std::unique_ptr<node> root_ptr = std::make_unique<node>();

        node *root = root_ptr.get();
        for (auto &p : values) {
            create_branch(root, p.first, 0, p.second);
        }
        return root_ptr;
    }

    int decode_command_with_huffman_tree(node *root) {
        if (root->val != -1) {
            return root->val;
        }
        bool next_bit = ch.get_bit();
        if (next_bit) {
            assert(root->right);
            return decode_command_with_huffman_tree(root->right);
        } else {
            assert(root->left);
            return decode_command_with_huffman_tree(root->left);
        }
    }

    void read_huffman_tree(node *root, int command_number, std::map<int, std::set<int>> &mapa) {
        int cur_length_symbol_pos = 0;
        int prev_command = -1;
        while (true) {
            if (cur_length_symbol_pos >= command_number) {
                break;
            }
            int command = decode_command_with_huffman_tree(root);
            assert(command >= 0 && command <= 18);
            int rep = -1;
            switch (command) {
                case 16:
                    assert(prev_command != -1);
                    rep = ch.get_n_bits_reversed_order(2);
                    rep += 3;
                    command = prev_command;
                case 17:
                    if (rep == -1) {
                        rep = ch.get_n_bits_reversed_order(3);
                        rep += 3;
                        command = 0;
                    }
                case 18:
                    if (rep == -1) {
                        rep = ch.get_n_bits_reversed_order(7);
                        rep += 11;
                        command = 0;
                    }
                    for (int i = 0; i < rep; ++i) {
                        assert(cur_length_symbol_pos < command_number);
                        mapa[command].insert(cur_length_symbol_pos++);
                    };
                    break;
                default:
                    mapa[command].insert(cur_length_symbol_pos++);
                    break;
            }
            prev_command = command;
        }
    }

    bool decode_dynamic_Huffman_codes() {
        int HLIT = ch.get_n_bits_reversed_order(5);
        int number_of_lengths_and_symbols = 257 + HLIT;
        int HDIST = ch.get_n_bits_reversed_order(5);
        int number_of_distances = 1 + HDIST;
        int HCLEN = ch.get_n_bits_reversed_order(4);//todo reverse?
        int number_of_lengths_for_first_huffman_tree = HCLEN + 4;
        std::vector<int> commands = {16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15};
        std::map<int, std::set<int>> tree_mapa;
        for (int i = 0; i < number_of_lengths_for_first_huffman_tree; ++i) {
            int k = ch.get_n_bits_reversed_order(3); // todo reversed
            tree_mapa[k].insert(commands[i]);
        }
        std::unique_ptr<node> root1 = create_Huffman_tree(tree_mapa);
        std::map<int, std::set<int>> length_symbols_mapa;
        read_huffman_tree(root1.get(), number_of_lengths_and_symbols, length_symbols_mapa);
        std::map<int, std::set<int>> length_distance_mapa;
        read_huffman_tree(root1.get(), number_of_distances, length_distance_mapa);
        std::unique_ptr<node> root_length = create_Huffman_tree(length_symbols_mapa);
        std::unique_ptr<node> root_distance = create_Huffman_tree(length_distance_mapa);
        auto get_command_func = [this, &root_length]() -> int {
            return decode_command_with_huffman_tree(root_length.get());
        };
        auto get_dist_func = [this, &root_distance]() -> int {
            return decode_command_with_huffman_tree(root_distance.get());
        };
        return decode_Huffman_codes(get_command_func, get_dist_func);
    }

    std::vector<decode_func_t> decode_funcions = {
            &zlib_decoder::decode_no_compression,
            &zlib_decoder::decode_fixed_Huffman_codes,
            &zlib_decoder::decode_dynamic_Huffman_codes
    };

    void write_to_buffer(uchar chh) {
        static std::ofstream dout("B:\\Projects\\GitProjects\\Graphics\\labs\\lab7\\debug.txt");
        static int kk = 0;
        static bool write_decoded_bytes = true;
        if (write_decoded_bytes) {
            if (cur_w == 0) {
                dout << "FILTER " << (int) chh << "\n";
                cur_w++;
            } else {
                if (byte_in_pixel == 1) {
                    dout << "color " << (int) chh << " pos w=" << cur_w++ << " pos h=" << cur_h << "\n";
                } else {
                    dout << "color " << (int) chh << " pos w[" << kk++ << "]=" << cur_w << " pos h=" << cur_h
                         << "\n";
                    if (kk == byte_in_pixel) {
                        ++cur_w;
                        kk = 0;
                    }
                }
                dout.flush();
                if (cur_w == max_w + 1) {
                    cur_w = 0;
                    ++cur_h;
                }
            }
        }
        if (output_buffer_pos >= max_buffer_size) {
            throw std::out_of_range("out of bounds for output image decoded buffer\n");
        }
        output_buffer[output_buffer_pos++] = chh;
    }

    chunk &ch;
    int CM = 0;
    int CINFO = 0;
    int FLEVEL = 0;
    int FDICT = 0;
    int FCHECK = 0;
    char *&output_buffer;
    size_t &output_buffer_pos;
    size_t max_buffer_size;
};

struct png_decoder {

    png_decoder(std::string file_name) : file_name(file_name) {};

    bool open_file() {
        assert(!file_name.empty());
        fin.open(input_file_name, std::ios::binary);
        if (!fin.is_open()) {
            std::cerr << "Failed to open input file \"" << input_file_name << "\"\n";
            return false;
        }
        return true;
    }

    bool read_signature() {
        uchar signature[] = {137, 80, 78, 71, 13, 10, 26, 10};
        for (int i = 0; i < 8; ++i) {
            uchar c = fin.get();
//            fin >> c;
            char k = char(0b11110000 & c);
            k >>= 4;
            char g = (char) (0b00001111 & c);
            if (c != signature[i]) {
                std::cerr << "Check for signature in input png file failed on position " << i << "\n";
                return false;
            }
        }
        return true;
    }

    bool read_IHDR() {
        chunk ch;
        ch.read(fin);
        if (ch.type != "IHDR") {
            std::cerr << "First chunk of input image is not IHDR chunk type. Input image incorrect.\n";
            return false;
        }
        if (ch.length != 13) {
            std::cerr << "Length of data in IHDR chunk is incorrect.\n";
            return false;
        }
        for (int i = 3; i >= 0; --i) {
            width = width | (ch.get_byte() << (8 * i));
        }
        max_w = width;
        for (int i = 3; i >= 0; --i) {
            height = height | (ch.get_byte() << (8 * i));
        }
        max_h = height;
        if (width * height == 0) {
            std::cerr << "Unsupported size (shouldn't be 0)\n";
            return false;
        }
        int depth = ch.get_byte();
        if (depth != 8) {
            std::cerr << "Unsupported bit depth = " << depth << ", should be 8\n";
            return false;
        }
        colorType = ch.get_byte();
        if (colorType != 0 && colorType != 2) {
            std::cerr << "Unsupported colorType = " << colorType << ", should be 0 or 2\n";
            return false;
        }
        byte_in_pixel = (colorType == 0) ? 1 : 3;
        int compressionMethod = ch.get_byte();
        if (compressionMethod != 0) {
            std::cerr << "Unsupported compressionMethod = " << compressionMethod << ", should be 0\n";
            return false;
        }
        filter_method = ch.get_byte();
        if (filter_method > 4 || filter_method < 0) {
            std::cerr << "Unsupported filter method = " << compressionMethod << ", should be in range [0...4]\n";
            return false;
        }
        int interlacedMethod = ch.get_byte();
        if (interlacedMethod) {
            std::cerr << "Program is not supporting interlacing\n";
            return false;
        }
        return true;
    }

    bool read_all_IDAT() {

        std::vector<chunk> v(1);
        v[0].read(fin);
        while (v[0].type != "IDAT") {
            chunk ch1;
            if (!ch1.read(fin)) {
                return false;
            }
            v[0].swap(ch1);
        }
        bool is_end = false;

        while (true) {
            std::cout << "Read IDAT " << v.back().length << "\n";
            v.emplace_back();
            v.back().read(fin);
            if (v.back().type != "IDAT") {
                if (v.back().type == "IEND"){
                    is_end = true;
                }
                v.pop_back();
                break;
            }
        }
        while (!is_end){
            chunk ch;
            ch.read(fin);
            if (ch.type == "IEND"){
                is_end = true;
            }
            if (ch.type == "IDAT") {
                std::cerr << "All chunks IDAT should appear consecutively with no other intervening chunks and"
                             " ends with IEND chunk\n";
                return false;
            }
        }
        char * input_buffer = nullptr;
        size_t input_buffer_size = 0;
        for (auto &i : v) {
            input_buffer_size += i.length;
        }
        try {
            input_buffer = new char[input_buffer_size];
        } catch (const std::bad_alloc &e) {
            std::cerr << "Allocation for chunk png archived data buffer failed: " << e.what() << '\n';
            return false;
        }
        size_t input_buffer_pos = 0;
        for (auto &cur_IDAT : v) {
            for (int i = 0; i < cur_IDAT.length; ++i) {
                input_buffer[input_buffer_pos++] = cur_IDAT.get_byte();
            }
        }
        chunk ch;
        ch.length = input_buffer_size;
        ch.data = input_buffer;
        pos = 0;
        buffer_size = width * height * byte_in_pixel + height;
        try {
            buffer = new char[buffer_size];
        } catch (const std::bad_alloc &e) {
            std::cerr << "Allocation for output data failed: " << e.what() << '\n';
            return false;
        }
        zlib_decoder decoder(ch, buffer, pos, buffer_size);
        if (!decoder.read_zlib_header()) {
            return false;
        }
        if (!decoder.decode_data()) {
            return false;
        }
        return true;
    }

    uchar get_char() {
        if (pos >= buffer_size) {
            throw std::out_of_range("out of range in output buffer for decoded data\n");
        }
        return (uchar) buffer[pos++]; // todo may delete?
    }

    static uchar none_filter(uchar x, uchar a, uchar b, uchar c) {
        return x;
    }

    static uchar sub_filter(uchar x, uchar a, uchar b, uchar c) {
        return x + a;
    }

    static uchar up_filter(uchar x, uchar a, uchar b, uchar c) {
        return x + b;
    }

    static uchar average_filter(uchar x, uchar a, uchar b, uchar c) {
        return x + (unsigned char) floor(((double) a + b) / 2.);
    }

    static uchar paeth_filter(uchar x, uchar a, uchar b, uchar c) {
        int p = a + b - c;
        int pa = abs(p - a);
        int pb = abs(p - b);
        int pc = abs(p - c);
        int Pr;
        if (pa <= pb and pa <= pc) {
            Pr = a;
        } else if (pb <= pc) {
            Pr = b;
        } else {
            Pr = c;
        }
        return x + Pr;
    }

    using filter_func_t = uchar (*)(uchar x, uchar a, uchar b, uchar c);

    std::vector<filter_func_t> filter_funcs = {
            &none_filter,
            &sub_filter,
            &up_filter,
            &average_filter,
            &paeth_filter
    };

    void process_data() {
        pos = 0;
        std::vector<uchar> prev_line(width *byte_in_pixel,
        0);
        for (int i = 0; i < height; ++i) {
            std::vector<uchar> cur_line(width *byte_in_pixel);
            uchar filter = get_char();
            assert(filter <= 4 && filter >= 0);
            std::vector<uchar> prev_element(byte_in_pixel, 0);
            for (int j = 0; j < width; ++j) {
                for (int k = 0; k < byte_in_pixel; ++k) {
                    uchar x = get_char();
                    uchar a = prev_element[k];
                    uchar b = prev_line[j * byte_in_pixel + k];
                    uchar c = (j == 0) ? 0 : prev_line[(j - 1) * byte_in_pixel + k];
                    uchar filtered_cur = filter_funcs[filter](x, a, b, c);
                    output.output << filtered_cur;
                    prev_element[k] = filtered_cur;
                    cur_line[j * byte_in_pixel + k] = filtered_cur;
                }
            }
            prev_line = cur_line;
        }
    }

    ~png_decoder() {
        delete[] buffer;
    }

    std::ifstream fin;
    std::string file_name;

    int width = 0;
    int height = 0;
    int colorType = 0;
    int filter_method = 0;
    char *buffer = nullptr;
    size_t buffer_size = 0;
    size_t pos = 0;
};


bool read_arguments(int argc, char *argv[]) {
    static std::string usage = "Usage : lab6.exe <name_of_input_file_1> <name_of_input_file_2> <name_of_output_files> <correlation_way>\n";
    if (argc != 3) {
        std::cerr << "Incorrect arguments number\n" << usage;
        return false;
    }
    input_file_name = argv[1];
    output_file_name = argv[2];
    return true;
}

void testing() {
    std::vector<std::string> res = {"0x15", "0x8d", "0x51", "0x0a", "0xc0", "0x20", "0x0c", "0x43", "0xff", "0x3d",
                                    "0x45", "0xae", "0x56", "0x67", "0xdd", "0x8a", "0x5d", "0x0b", "0xd5", "0x21",
                                    "0xde", "0x7e", "0x0a", "0xf9", "0x08", "0x21", "0x2f", "0xc9", "0x4a", "0x57",
                                    "0xcb", "0x12", "0x05", "0x5d", "0xec", "0xde", "0x82", "0x18", "0xc6", "0xc3",
                                    "0x28", "0x4c", "0x05", "0x5e", "0x61", "0x72", "0x3f", "0x23", "0x0d", "0x6a",
                                    "0x7c", "0xe2", "0xce", "0xc8", "0xe1", "0x8d", "0x0d", "0x73", "0x77", "0x3b",
                                    "0xc8", "0x0a", "0x94", "0x29", "0x36", "0xe3", "0xa8", "0xba", "0x12", "0xa9",
                                    "0x62", "0xf9", "0x17", "0x50", "0xa9", "0x9c", "0xb6", "0xc3", "0xe4", "0x60",
                                    "0xb8", "0xe9", "0xc2", "0x24", "0x19", "0xe7", "0xa1", "0x7a", "0xec", "0x2d",
                                    "0xe9", "0x78", "0xfd", "0x65", "0x1b", "0x07", "0xa5", "0x90", "0xce", "0xe9",
                                    "0x07"};
    std::vector<std::string> res2 = {"05", "A0", "61", "0B", "00", "00", "18", "86", "2E", "68", "8F", "EA", "D9", "3D",
                                     "AE", "C7", "B4", "9B", "10", "E8", "FF", "40", "21", "07", "14", "56", "5A",
                                     "DA"};

    std::vector<std::string> res3 = {"63", "F8", "3F", "93", "E1", "3F", "03", "C3", "CC", "FF", "20", "1A", "C8", "00",
                                     "22", "24", "0E", "58", "12", "85", "33", "D3", "F8", "3F", "03", "32", "07", "44",
                                     "03", "00", "AA", "05", "23", "77"};
    int k = 0;
    char *buffer = new char[res.size()];
    for (int i = 0; i < res.size(); ++i) {
        std::string s = res[i];
        std::stringstream ss;
        ss << std::hex << s;
        unsigned n;
        ss >> n;
        buffer[i] = n;
    }
    chunk ch;
    ch.data = buffer;
    ch.length = res.size();
    size_t max_output_buffer_size = res.size() * 100;
    char *output_buffer = new char[max_output_buffer_size];
    size_t pos = 0;
    zlib_decoder decoder(ch, output_buffer, pos, max_output_buffer_size);
    decoder.decode_data();
    int gg = 0;
    std::cout << "\n" << output_buffer;

    /* for (int i = 0; i < 4; i++){
         std::cout << (int) ((uchar) output_buffer[gg++] )<< "\n";
         for (int j = 0; j < 4; ++j){
             for (int t = 0; t < 4; ++t){
                 std::cout << (int) ((uchar) output_buffer[gg++] )<< " ";
             }
             std::cout <<"\n";
         }
     }*/


//    std::cout << "HERE:\n" << output_buffer;

}

void gen_image() {
    file_worker wo;
    wo.output_file_name = "B:\\Projects\\GitProjects\\Graphics\\pictures\\source_images\\small_rectangle.ppm";
    int w = 40;
    int h = 20;
    wo.open_output_file();
    wo.set_header_for_output_file(false, w, h, 255);
    wo.write_header();

    for (int i = 0; i < h; ++i) {
        for (int j = 0; j < w; ++j) {
            wo.output << uchar(j * 30) << uchar(i * 50) << uchar((j + i) * 10);
        }
    }
}

int main(int argc, char *argv[]) {
//    testing();
//    gen_image();
//    return 0;
    bool testing = true;
    std::vector<std::string> names = {
            "seeds_no_alpha",
//            "sm_rec",
//            "site_ex",
    };
    int i = names.size() - 1;
    if (testing) {
        std::string name = names[i];
        int argct = 3;
        char **argvt = new char *[argct];
        argvt[0] = "lab3.exe";
        std::string name_file = "B:\\Projects\\GitProjects\\Graphics\\pictures\\source_images\\" + name + ".png";
        argvt[1] = const_cast<char *>(name_file.c_str());
        std::string name_file_out =
                "B:\\Projects\\GitProjects\\Graphics\\pictures\\output_pictures\\" + name + ".ppm";
        argvt[2] = const_cast<char *>(name_file_out.c_str());
        if (!read_arguments(argct, argvt)) {
            return 1;
        }
    } else {
        if (!read_arguments(argc, argv)) {
            return 1;
        }
    }
    png_decoder decoder(input_file_name);
    if (!decoder.open_file()) {
        return 1;
    }
    if (!decoder.read_signature()) {
        return 1;
    }
    if (!decoder.read_IHDR()) {
        return 1;
    }

    if (!decoder.read_all_IDAT()) {
        return 1;
    }
    output.output_file_name = output_file_name;
    if (!output.open_output_file()) {
        return 1;
    }
    output.set_header_for_output_file(decoder.colorType == 0, decoder.width, decoder.height, 255);
    output.write_header();
    decoder.process_data();
    return 0;
}