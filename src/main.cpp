#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>

int32_t ABS(int32_t a) {
	return (a > 0) ? a : -1 * a;
}

uint32_t LUT_base_lengths[] = {
	3, 4, 5, 6, 7, 8, 9, 10, // 257 - 264
	11, 13, 15, 17,          // 265 - 268
	19, 23, 27, 31,          // 269 - 273
	35, 43, 51, 59,          // 274 - 276
	67, 83, 99, 115,         // 278 - 280
	131, 173, 195, 227,      // 281 - 284
	258						 // 285
};

uint8_t LUT_base_lengths_extra_bit[] = {
	0, 0, 0, 0, 0, 0, 0, 0, // 257 - 264
	1, 1, 1, 1,             // 265 - 268
	2, 2, 2, 2,             // 269 - 273
	3, 3, 3, 3,             // 274 - 276
	4, 4, 4, 4,             // 278 - 280
	5, 5, 5, 5,             // 281 - 284
	0						// 285
};

uint32_t LUT_dist_bases[] = {
    /*0*/ 1, 2, 3, 4,    //0-3
    /*1*/ 5, 7,          //4-5
    /*2*/ 9, 13,         //6-7
    /*3*/ 17, 25,        //8-9
    /*4*/ 33, 49,        //10-11
    /*5*/ 65, 97,        //12-13
    /*6*/ 129, 193,      //14-15
    /*7*/ 257, 385,      //16-17
    /*8*/ 513, 769,      //18-19
    /*9*/ 1025, 1537,    //20-21
    /*10*/ 2049, 3073,   //22-23
    /*11*/ 4097, 6145,   //24-25
    /*12*/ 8193, 12289,  //26-27
    /*13*/ 16385, 24577, //28-29
             0   , 0      //30-31, error, shouldn't occur
};

uint32_t LUT_dist_extra_bits[] = {
    /*0*/ 0, 0, 0, 0, //0-3
    /*1*/ 1, 1,       //4-5
    /*2*/ 2, 2,       //6-7
    /*3*/ 3, 3,       //8-9
    /*4*/ 4, 4,       //10-11
    /*5*/ 5, 5,       //12-13
    /*6*/ 6, 6,       //14-15
    /*7*/ 7, 7,       //16-17
    /*8*/ 8, 8,       //18-19
    /*9*/ 9, 9,       //20-21
    /*10*/ 10, 10,    //22-23
    /*11*/ 11, 11,    //24-25
    /*12*/ 12, 12,    //26-27
    /*13*/ 13, 13,     //28-29
            0 , 0      //30-31 error, they shouldn't occur
};

uint8_t LUT_default_lengths[] = {
   8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8, 8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
   8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8, 8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
   8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8, 8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
   8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8, 8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
   8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8, 9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,
   9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9, 9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,
   9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9, 9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,
   9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9, 9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,
   7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7, 7,7,7,7,7,7,7,7,8,8,8,8,8,8,8,8
};

uint8_t LUT_default_distance[] = {
   5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5
};

uint32_t EndianSwap(uint32_t a) {
	uint32_t b;
	b = (a << 24 |
		(a << 8) & 0x00FF0000 |
		(a >> 8) & 0x0000FF00 |
		 a >> 24 );
	return b;
};

struct mem_buffer {
	uint32_t size;
	void* buffer;
};

struct seek_buffer {
	void* buffer;
	uint32_t size;
	uint32_t bytes_read;
	uint32_t bits_read;
};

struct sliding_buffer {
	uint8_t* buffer;
	uint32_t total_size;
	uint32_t actual_size;
};

bool CheckBounds(seek_buffer* sb, uint32_t no_of_bits) {
	if((sb->bytes_read + (sb->bits_read + no_of_bits) / 8) > sb->size) {
		//printf("ERROR: Buffer Overflow \n");
		return false;
	}
	return true;
}

void FlushByte(seek_buffer* sb) {
	if(sb->bits_read) {
		sb->bytes_read++;
		sb->bits_read = 0;
	}
}

void SeekBuffer(seek_buffer* sb, uint32_t n) {
	if(CheckBounds(sb, n)) {
		sb->bits_read += n;
		uint32_t a = sb->bits_read % 8;
		uint32_t b = sb->bits_read / 8;
		sb->bytes_read += b;
		sb->bits_read = a;
	}
}

uint32_t PeekBits(seek_buffer* sb, uint32_t n) {
	uint32_t result = 0;
	if(CheckBounds(sb, n)) {
		uint32_t mask = ((1 << n) - 1) << sb->bits_read;
		result = *((uint32_t*)((uint8_t*)sb->buffer + sb->bytes_read)) & mask;
		result = result >> sb->bits_read;
	}
	return result;
}

uint32_t ExtractBits(seek_buffer* sb, uint32_t n) {
	uint32_t result = PeekBits(sb, n);
	SeekBuffer(sb, n);
	return result;
}

uint32_t ReverseBits(uint32_t bits, uint32_t n) {
	uint32_t result = 0;
	for(uint32_t i=0; i<n; i++) {
		result <<= 1;
		uint32_t bit = bits & (1<<i);
		result |= (bit>0) ? 1 : 0;
	}
	return result;
}

struct zhuffman {
	uint32_t* assigned_codes;
	uint8_t* codelengths;
	uint32_t n;
};

struct zlib_stream {
	void* buffer;
	uint32_t size;

	zlib_stream* next;
};

mem_buffer ExtractBufferFromLinkedList(zlib_stream* list) { 
	mem_buffer result = {};
	uint32_t total_size = 0;
	zlib_stream* temp = list;
	for(;;) {
		total_size += temp->size;
		if(temp->next) temp = temp->next;
		else break;
	}
	result.size = total_size;
	temp = list;
	void* new_buffer = malloc(total_size);
	void* seek = new_buffer;
	while(total_size != 0 ) {
		memcpy(seek , temp->buffer, temp->size);
		seek = (void*)((uint8_t*)seek + temp->size);
		total_size -= temp->size;
		temp = temp->next;
	}
	result.buffer = new_buffer;
	return result;
}

void zlib_stream_add_end(zlib_stream* list, void* buffer, uint32_t size) {
	if(!list->buffer) {
		list->buffer = malloc(size);
		list->size = size;
		memcpy(list->buffer, buffer, size);
		list->next = NULL;
	} else {
		zlib_stream* new_stream = (zlib_stream*) malloc(sizeof(zlib_stream));
		new_stream->buffer = malloc(size);
		new_stream->size = size;
		new_stream->next = NULL;
		memcpy(new_stream->buffer, buffer, size);

		zlib_stream* temp = list;
		while(temp->next != NULL) {
			temp = temp->next;
		}
		temp->next = new_stream;
	}
}

#pragma pack(push, 1)
int PNG_FILE_SIGNATURE[8] = { 137, 80, 78, 71, 13, 10, 26, 10 };

struct png_chunk_header {
	uint32_t length;
	uint32_t type;
};

struct png_chunk_footer {
	uint32_t crc;
};

struct png_ihdr {
	uint32_t width;
	uint32_t height;
	uint8_t bit_depth;
	uint8_t color_type;
	uint8_t compression_method;
	uint8_t filter_method;
	uint8_t interlace_method;
};

struct zlib_stream_header{
	uint8_t compression_flags;
	uint8_t additional_flags;
};

struct zlib_stream_footer {
	uint32_t check_value;
};

#pragma pack(pop)

enum filter_types {
	NO_FILTER,
	SUB_FILTER,
	UP_FILTER,
	AVG_FILTER,
	PAETH_FILTER
};

#pragma pack(push, 1)
struct bitmap_file_header {
	uint16_t type;
	uint32_t size;
	uint16_t reserved1;
	uint16_t reserved2;
	uint32_t offset;
};

struct bitmap_info_header {
	uint32_t size;
	int32_t width;
	int32_t height;
	uint16_t planes;
	uint16_t bitcount;
	uint32_t compression;
	uint32_t image_size;
	int32_t hor_ppm;
	int32_t ver_ppm;
	uint32_t colors_in_palette;
	uint32_t imp_colors;
};

#pragma pack(pop)

struct final_image {
	uint8_t* pixels;
	uint32_t size;
	uint32_t width;
	uint32_t height;
};

mem_buffer ReadEntireFile(char* Filename) {
	mem_buffer Result = {};
	
	FILE* File = fopen(Filename, "rb");
	if(File) {
		fseek(File, 0, SEEK_END);
		Result.size = ftell(File);
		fseek(File, 0, SEEK_SET);
		Result.buffer = malloc(Result.size);
		fread(Result.buffer, Result.size, 1, File);
		fclose(File);
		//printf("\nINFO: FileName: %s \n", Filename);
		//printf("INFO: FileSize: %d \n", Result.size);
	} else {
		//printf("FAILED: Cannot open file %s \n", Filename);
	}

	return Result;
}

uint32_t* BuildHuffman(uint8_t* codelengths, uint32_t n) {
	uint8_t max_bit_length = 0;
	for(uint32_t i=0; i<n; i++) 
		if(max_bit_length < codelengths[i]) max_bit_length = codelengths[i]; 

	uint32_t *code_count = (uint32_t*) calloc(max_bit_length+1, sizeof(uint32_t));
	uint32_t *first_codes = (uint32_t*) calloc(max_bit_length+1, sizeof(uint32_t));
	uint32_t *assigned_codes = (uint32_t*) calloc(n, sizeof(uint32_t));

	for(uint32_t i=0; i<n; i++) 
		code_count[codelengths[i]]++;
	code_count[0] = 0;

	uint32_t code = 0;
	for(uint32_t i=1; i<=max_bit_length; i++) {
		code = (code + code_count[i-1]) << 1;
		if(code_count[i] > 0) first_codes[i] = code;
	}

	for(uint32_t i=0; i<n; i++) 
		if(codelengths[i]) assigned_codes[i] = first_codes[codelengths[i]]++;

	return assigned_codes;
}

uint32_t DecodeHuffman(seek_buffer* sb, uint32_t* assigned_codes, uint8_t* codelengths, uint32_t n) {
	for(uint32_t i=0; i<n; i++) {
		if(codelengths[i] == 0) continue;
		uint32_t code = PeekBits(sb, codelengths[i]);
		code = ReverseBits(code, codelengths[i]);
		if(assigned_codes[i] == code) {
			SeekBuffer(sb, codelengths[i]); 
			return i;
		}
	}
	return 0;
}

void ParseHuffman(sliding_buffer* mem, seek_buffer* sb, zhuffman* lit_len, zhuffman* dist) {
	uint32_t data_index = 0;
	for(;;) {
		uint32_t decoded_value = DecodeHuffman(sb, lit_len->assigned_codes, lit_len->codelengths, lit_len->n);
		if(decoded_value == 256) break; 
		if(decoded_value < 256) {
			mem->buffer[mem->actual_size + data_index++] = decoded_value;
			continue;
		}
		if(decoded_value > 256 && decoded_value < 286) {
			uint32_t base_index = decoded_value - 257;
			uint32_t duplicate_length = LUT_base_lengths[base_index] +  ExtractBits(sb, LUT_base_lengths_extra_bit[base_index]);

			uint32_t distance_index = DecodeHuffman(sb, dist->assigned_codes, dist->codelengths, dist->n);
			uint32_t distance_length = LUT_dist_bases[distance_index] + ExtractBits(sb, LUT_dist_extra_bits[distance_index]);

			uint32_t back_pointer = mem->actual_size + data_index - distance_length;
			while(duplicate_length--) {
				mem->buffer[mem->actual_size + data_index++] = mem->buffer[back_pointer++];
			}
		}
	}
	mem->actual_size += data_index;
}

uint8_t PaethPredict(int32_t a, int32_t b, int32_t c) {
	int32_t p = a + b + c;
	int32_t pa = ABS(p-a); 
	int32_t pb = ABS(p-b); 
	int32_t pc = ABS(p-c); 

	if(pa<=pb && pa<=pc) return a;
	if(pb<=pc) return b;
	return c;
}

mem_buffer DecompressImage(zlib_stream* stream) {
	mem_buffer result = {};
	mem_buffer local_data = ExtractBufferFromLinkedList(stream);
	seek_buffer sb = {};
	sb.buffer = local_data.buffer;
	sb.size = local_data.size;

	zlib_stream_header* header = (zlib_stream_header*)sb.buffer;
	sb.bytes_read += sizeof(zlib_stream_header);
	uint8_t CM     = header->compression_flags & 0x0F;
	uint8_t CINFO  = header->compression_flags >> 4;
	uint8_t FCHECK = header->additional_flags & 0x1F;
	uint8_t FDICT  = (header->additional_flags >> 5) & 0x01;
	uint8_t FLEVEL = header->additional_flags >> 6;

	//printf("INFO: Compression Method - %u\n", CM);
	//printf("INFO: Compression Info   - %u\n", CINFO);
	//printf("INFO: FCHECK             - %u\n", FCHECK);
	//printf("INFO: FDICT              - %u\n", FDICT);
	//printf("INFO: FLEVEL             - %u\n", FLEVEL);

	if(CM != 8) {
		//printf("FAILED: Invalid compression method\n");
		return result;
	} else if(FDICT == 1) {
		//printf("FAILED: Invalid compression method - Preset Dictionary \n");
		return result;
	}

	uint32_t data_read = 0;
	CM = header->compression_flags & 0x0F;
	sliding_buffer decompressed_data = {};
	decompressed_data.total_size = 1024*1024*4;
	decompressed_data.buffer = (uint8_t*) malloc(decompressed_data.total_size);

	zhuffman lit_len;
	zhuffman dist;

	uint32_t BFINAL = 0;
	while(BFINAL == 0) {
		BFINAL = ExtractBits(&sb, 1);
		uint32_t BTYPE = ExtractBits(&sb, 2);

		//printf("INFO: BFINAL - %u\n", BFINAL);
		//printf("INFO: BTYPE  - %u\n", BTYPE);

		if(BTYPE == 0) {
			//printf("\nLOG: Compression Method: None\n");
			FlushByte(&sb);
			uint16_t LEN = ExtractBits(&sb, 16);
			uint16_t NLEN = ExtractBits(&sb, 16);

			//printf("INFO: LEN  - %u\n", LEN);
			//printf("INFO: NLEN  - %u\n", NLEN);
			//printf("INFO: ~LEN  - %u\n", ~LEN);
			//printf("INFO: ~NLEN  - %u\n", ~NLEN);

			if((uint16_t)LEN != ~(uint16_t)NLEN) {
				//printf("FAILED: Corrupted data, NLEN is not ones complement of LEN");
				return result;
			}

			return result;

		} else if(BTYPE == 1) {
			//printf("\nLOG: Compression Method: Fixed Huffman code\n");
			lit_len.codelengths = LUT_default_lengths;
			lit_len.n = 288;
			lit_len.assigned_codes = BuildHuffman(LUT_default_lengths, 288);

			dist.codelengths = LUT_default_distance;
			dist.n = 32;
			dist.assigned_codes = BuildHuffman(LUT_default_distance, 32);

		} else if (BTYPE == 2) {
			//printf("\nLOG: Compression Method: Dynamic Huffman code\n");
			uint32_t HLIT = ExtractBits(&sb, 5);
			uint32_t HDIST = ExtractBits(&sb, 5);
			uint32_t HCLEN = ExtractBits(&sb, 4);

			HLIT += 257;
			HDIST += 1;
			HCLEN += 4;

			//printf("\nLOG: HLIT: %u\n", HLIT);
			//printf("LOG: HDIST: %u\n", HDIST);
			//printf("LOG: HCLEN: %u\n", HCLEN);

			uint32_t HCLEN_swizzle[] = { 16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15 };
			uint8_t HCLEN_table[19] = {};
			for(uint32_t i = 0; i<HCLEN; i++) {
				HCLEN_table[HCLEN_swizzle[i]] = ExtractBits(&sb, 3);
			}

			uint32_t* tree_of_two_trees = BuildHuffman(HCLEN_table, 19); 
			uint8_t* two_trees_codelengths = (uint8_t*) malloc(HLIT+HDIST);

			uint32_t code_index = 0;
			while(code_index < (HDIST+HLIT)) {
				uint32_t decoded_value = DecodeHuffman(&sb, tree_of_two_trees, HCLEN_table, 19);
				if(decoded_value < 16) {
					two_trees_codelengths[code_index++] = decoded_value;
					continue;
				}
				uint32_t repeat_count = 0;
				uint8_t code_length = 0;

				switch(decoded_value) {
					case 16:
						repeat_count = ExtractBits(&sb, 2) + 3;
						code_length = two_trees_codelengths[code_index - 1];
						break;
					case 17:
						repeat_count = ExtractBits(&sb, 3) + 3;
						break;
					case 18:
						repeat_count = ExtractBits(&sb, 7) + 11;
						break;
				}
				memset(two_trees_codelengths + code_index, code_length, repeat_count); 
				code_index += repeat_count;
			}

			lit_len.assigned_codes =  BuildHuffman(two_trees_codelengths, HLIT);
			lit_len.codelengths = two_trees_codelengths;
			lit_len.n = HLIT;

			dist.assigned_codes = BuildHuffman(two_trees_codelengths + HLIT, HDIST);
			dist.codelengths = two_trees_codelengths + HLIT;
			dist.n = HDIST;
		}
		ParseHuffman(&decompressed_data, &sb, &lit_len, &dist); 
	}
	result.size = decompressed_data.actual_size;
	result.buffer = decompressed_data.buffer;
	return result;
}

mem_buffer DefilterImage(uint8_t* data, png_ihdr* ihdr) {
	mem_buffer result = {};
	uint32_t x = ihdr->width;
	uint32_t y = ihdr->height;
	uint8_t bit_depth = ihdr->bit_depth;
	uint8_t bytes_per_pixel = 4;

	uint8_t* row = data;
	uint32_t stride = x * bytes_per_pixel;
	uint8_t* image = (uint8_t*) malloc(x*y*bytes_per_pixel);
	uint8_t* working = image;

	for(uint32_t i=0; i<y; i++) {
		working = image + i*stride;
		uint8_t filter = *row++;
		uint8_t* prev_row = working - stride;

		switch(filter) {
			case NO_FILTER:
				for(uint32_t j=0; j<x; j++) {
					for(uint32_t k=0; k<bytes_per_pixel; k++) {
						working[j*bytes_per_pixel + k] = row[j*bytes_per_pixel +k];
					}
				} break;

			case SUB_FILTER:
				for(uint32_t j=0; j<x; j++) {
					for(uint32_t k=0; k<bytes_per_pixel; k++) {
						uint8_t a =0;
						if(j) a = working[(j-1)*bytes_per_pixel + k];
						uint8_t value = row[j*bytes_per_pixel + k] + a;
						working[j*bytes_per_pixel + k] = value;
					}
				} break;

			case AVG_FILTER:
				for(uint32_t j=0; j<x; j++) {
					for(uint32_t k=0; k<bytes_per_pixel; k++) {
						uint8_t a = 0;
						uint8_t b = prev_row[j+bytes_per_pixel + k];
						if(j) a = working[(j-1)*bytes_per_pixel + k];
						uint8_t value = row[j*bytes_per_pixel + k] + ((a+b) >> 1);
						working[j*bytes_per_pixel + k] = value;
					}
				} break;

			case PAETH_FILTER:
				for(uint32_t j=0; j<x; j++) {
					for(uint32_t k=0; k<bytes_per_pixel; k++) {
						uint8_t a = 0;
						uint8_t b = prev_row[j*bytes_per_pixel + k];
						uint8_t c = 0;
						if(j) {
							a = working[(j-1)*bytes_per_pixel + k];
							a = prev_row[(j-1)*bytes_per_pixel + k];
						}
						uint8_t value = row[j] + PaethPredict((int32_t) a, (int32_t) b, (int32_t) c);
						working[j*bytes_per_pixel + k] = value;
					}
				} break;
		}
		row += stride;
	}
	result.buffer = image;
	result.size = x*y*bytes_per_pixel;
	return result;
}

bool ParsePNG(final_image* image, mem_buffer* File) {
	uint8_t* buffer =(uint8_t*) File->buffer;
	uint32_t offset = 0;
	bool read_ihdr = false;
	png_ihdr* ihdr;
	uint32_t count = 0;
	zlib_stream zlib_compressed_data = {0};

	for(offset = 0; offset < 8; offset ++) {
		if(buffer[offset] != PNG_FILE_SIGNATURE[offset ]) {
			//printf("\nFAILED: Missing/bad PNG file signature \n");
			return false;
		}
	}
	//printf("\nPASSED: Read PNG File signature\n");
	
	for(;;) {
		png_chunk_header* header = (png_chunk_header*)(buffer + offset);
		uint32_t chunk_length = EndianSwap(header->length);
		uint32_t chunk_type   = EndianSwap(header->type);
		void* chunk_data   = (void*)(header + 1);

		png_chunk_footer* footer = (png_chunk_footer*)((uint8_t*)chunk_data + chunk_length);
		uint32_t chunk_crc    = EndianSwap(footer->crc);
		offset += chunk_length + sizeof(png_chunk_header) + sizeof(png_chunk_footer);
		//printf("\nINFO: Chunk type   - %c%c%c%c \n", 
										//((uint8_t*)(&chunk_type))[3],
										//((uint8_t*)(&chunk_type))[2],
										//((uint8_t*)(&chunk_type))[1],
										//((uint8_t*)(&chunk_type))[0] );
		//printf("\nINFO: Chunk length - %u \n", chunk_length);
		if(!read_ihdr) {
			if(chunk_type != 'IHDR') {
				//printf("\n FAILED: Missing IHDR Chunk \n");
				return false;
			}
			read_ihdr = true;
			//printf("\nLOG: Reading IHDR Chunk \n");
			ihdr = (png_ihdr*)chunk_data;
			ihdr->width = EndianSwap(ihdr->width);
			ihdr->height = EndianSwap(ihdr->height);

			//printf("INFO: Width              - %u \n",ihdr->width);
			//printf("INFO: Height             - %u \n",ihdr->height);
			//printf("INFO: Bit Depth          - %u \n",ihdr->bit_depth);
			//printf("INFO: Color Type         - %u \n",ihdr->color_type);
			//printf("INFO: Compression Method - %u \n",ihdr->compression_method);
			//printf("INFO: Filter Method      - %u \n",ihdr->filter_method);
			//printf("INFO: Interlace Method   - %u \n",ihdr->interlace_method);

			if(ihdr->compression_method != 0) {
				//printf("\n FAILED: Invalid compression method \n");
				return false;
			}
			if(ihdr->filter_method != 0) {
				//printf("\n FAILED: Invalid filter method \n");
				return false;
			}
			if(ihdr->interlace_method != 0) {
				//printf("\n FAILED: Unsupported interlace method \n");
				return false;
			}
			if(ihdr->color_type != 6) {
				//printf("\n FAILED: Unsupported color type\n");
				return false;
			}

		} else if(chunk_type =='sRGB') {
			//printf("\nWARNING: Encountered sRGB chunk, ignoring ancillary chunk\n");

		} else if(chunk_type == 'IDAT') {
			zlib_stream_add_end(&zlib_compressed_data, chunk_data, chunk_length);

		} else if(chunk_type == 'IEND') {
			//printf("\nLOG: Hit IEND Chunk \n");
			mem_buffer decompressed_data = DecompressImage(&zlib_compressed_data);
			mem_buffer reconstructed_image = DefilterImage((uint8_t*) decompressed_data.buffer, ihdr);

			image->pixels = (uint8_t*) decompressed_data.buffer;
			image->size = decompressed_data.size;
			image->width = ihdr->width;
			image->height = ihdr->height;
			
			return true;

		} else {
			//printf("\nFAILED: Unsupported/invalid Chunk \n");
		}
	}
}

void ConstructBMPPixelArray(final_image* image) {
	uint32_t padded_row_size = ((24*image->width + 31)/32) * 4;
	uint8_t* buffer = (uint8_t*) calloc(1, padded_row_size * image->height);

	int32_t row_index_bmp = image->height - 1;

	uint32_t row_index_png = 0;
	uint32_t row_length_png = 4*image->width;

	while(row_index_bmp >= 0) {
		uint32_t pixel_index = 0;
		while(pixel_index < image->width) {
			uint32_t color_index = 0;
			while(color_index < 3) {
				buffer[row_index_bmp*padded_row_size + pixel_index*3 + color_index] = 
					image->pixels[row_index_png*row_length_png + pixel_index*4 + color_index];
				color_index++;
			}
			pixel_index++;
		}
		row_index_bmp--;
		row_index_png++;
	}
	
	image->pixels = buffer;
	image->size = padded_row_size * image->height;
}

void ConstructBMP(final_image* image) {
	ConstructBMPPixelArray(image);

	bitmap_file_header file_header = {};
	file_header.type = 'MB';
	file_header.size = sizeof(bitmap_file_header) + sizeof(bitmap_info_header) + image->size;
	file_header.reserved1 = 0;
	file_header.reserved2 = 0;
	file_header.offset = sizeof(bitmap_file_header) + sizeof(bitmap_info_header);
	
	bitmap_info_header info_header = {};
	info_header.size = 40;
	info_header.width = image->width;
	info_header.height = image->height;
	info_header.planes = 1;
	info_header.bitcount = 24;
	info_header.compression = 0;
	info_header.image_size = 0;
	info_header.hor_ppm = 0;
	info_header.ver_ppm = 0;
	info_header.colors_in_palette = 0;
	info_header.imp_colors = 0;

	FILE* file;
	file = fopen("output.bmp", "w");
	fwrite(&file_header, sizeof(bitmap_file_header), 1, file);
	fwrite(&info_header, sizeof(bitmap_info_header), 1, file);
	fwrite(image->pixels, image->size, 1, file);
}

void ConstructPPM(final_image* image) {
	std::cout << "P3\n" << image->width << ' ' << image->height << "\n255\n";
	uint32_t pixelcount = image->width * image->height;
	for(uint32_t i =0; i<pixelcount; i++) {
		std::cout << (int)image->pixels[i*4] << ' ' << (int)image->pixels[i*4 + 1] << ' ' << (int)image->pixels[i*4 + 2] << '\n';
	}
}

int main(int argc, char** argv) {

	mem_buffer File = ReadEntireFile("test.png");

	final_image img = {};
	ParsePNG(&img, &File);
	//ConstructBMP(&img);
	ConstructPPM(&img);

	return 0;
}
