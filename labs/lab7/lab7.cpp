#include <iostream>
#include <fstream>
#include <utility>
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

std::string input_file_name;
std::string output_file_name;
using uchar = unsigned char;

int byte_in_pixel;

struct file_worker {
    explicit file_worker() = default;

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

    void write_header() {
        assert(output.is_open());
        output << (is_P5 ? "P5" : "P6") << '\n';
        output << std::to_string(width) + " " + std::to_string(height) << '\n';
        output << std::to_string(max_color_size) << '\n';
    }

    std::string output_file_name;

    char *input_buffer = nullptr;
    std::ifstream input;
    std::ofstream output;

    size_t width;
    size_t height;
    size_t max_color_size;

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
            cur_byte = (uchar) get_byte();
            pos_in_byte = 0;
        }
        return cur_byte & (1 << (pos_in_byte++));
    }

    uchar get_byte() {
        assert(pos_in_byte == 8);
        if (pos >= length) {
            throw std::out_of_range("no readable elements in chunk");
        }
        return data[pos++];
    }

    int get_n_bits_reversed_order(int n) {
        int res = 0;
        for (int i = 0; i < n; i++) {
            bool ch = get_bit();
            res |= (ch << i);
        }
        return res;
    }

    int get_n_bits_direct_order(int n) {
        int res = 0;
        for (int i = 0; i < n; i++) {
            bool chh = get_bit();
            res <<= 1;
            res |= chh;
        }
        return res;
    }


    uint32_t length = 0;
    std::string type = "1234";
    char *data = nullptr;
    uchar CRC[4];
    size_t pos = 0;
    uchar cur_byte;
    int pos_in_byte = 8;
};

struct zlib_decoder {

    zlib_decoder(chunk &ch, char *&outputBuffer, size_t &outputBufferPos, size_t max_buffer_size)
            : ch(ch), output_buffer(outputBuffer), output_buffer_pos(outputBufferPos),
              max_buffer_size(max_buffer_size) {}

    bool read_zlib_header() {
        uchar CMF = ch.get_byte();
        CM = 0b00001111 & CMF;
        if (CM != 8) {
            std::cerr << "CM in zlib header != 8\n";
            return false;
        }
        uchar FLG = ch.get_byte();
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
        uchar b1 = ch.get_byte();
        uchar b2 = ch.get_byte();
        int length = (b1 << 8) | b2;
        uchar b1c = ch.get_byte();
        uchar b2c = ch.get_byte();
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
                if (add_bits != 0) {
                    shift += ch.get_n_bits_reversed_order(add_bits);;
                }
                size_t initial_buffer_pos = output_buffer_pos;
                size_t shifted_position = initial_buffer_pos - shift;
                for (int j = 0; j < length; j++) {
                    if (shifted_position == initial_buffer_pos) {
                        shifted_position = initial_buffer_pos - shift;
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
        int HCLEN = ch.get_n_bits_reversed_order(4);
        int number_of_lengths_for_first_huffman_tree = HCLEN + 4;
        std::vector<int> commands = {16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15};
        std::map<int, std::set<int>> tree_mapa;
        for (int i = 0; i < number_of_lengths_for_first_huffman_tree; ++i) {
            int k = ch.get_n_bits_reversed_order(3);
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
        if (output_buffer_pos >= max_buffer_size) {
            throw std::out_of_range("out of bounds for output image decoded buffer\n");
        }
        output_buffer[output_buffer_pos++] = (char) chh;
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

    png_decoder(std::string file_name) : file_name(std::move(file_name)) {};

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
        for (int i = 3; i >= 0; --i) {
            height = height | (ch.get_byte() << (8 * i));
        }
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
            v.emplace_back();
            v.back().read(fin);
            if (v.back().type != "IDAT") {
                if (v.back().type == "IEND") {
                    is_end = true;
                }
                v.pop_back();
                break;
            }
        }
        while (!is_end) {
            chunk ch;
            ch.read(fin);
            if (ch.type == "IEND") {
                is_end = true;
            }
            if (ch.type == "IDAT") {
                std::cerr << "All chunks IDAT should appear consecutively with no other intervening chunks and"
                             " ends with IEND chunk\n";
                return false;
            }
        }
        char *input_buffer = nullptr;
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
                input_buffer[input_buffer_pos++] = (char) cur_IDAT.get_byte();
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
        return (uchar) buffer[pos++];
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
        return x + (uchar) floor(((double) a + b) / 2.);
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


int main(int argc, char *argv[]) {
    if (!read_arguments(argc, argv)) {
        return 1;
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