#include <png.h>
#include <stdio.h>
#include <stdlib.h>

int PlanePNG(char* folder, int timestep, int ncol, int nrow, int size, unsigned char* data) {
	int i,j;
	char name[512];
	FILE *fp;
	static png_structp png_ptr;
	static png_infop info_ptr;
	static png_bytep ptr;
	
	sprintf(name,"%s/role_%d.png", folder,timestep);
//	printf("name: %s\n", name);	
	
	fp = fopen(name,"wb");
	png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, (png_voidp)NULL, (png_error_ptr)NULL, (png_error_ptr)NULL);
	info_ptr = png_create_info_struct (png_ptr);
	png_init_io(png_ptr, fp);
	png_set_IHDR(png_ptr, info_ptr, ncol, nrow, 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
	png_write_info(png_ptr,info_ptr);
	
	for (i=0; i < nrow; i++) {
		ptr = data + i*3*ncol;
		png_write_rows(png_ptr, &ptr, 1);
	}
	
	png_write_end(png_ptr, info_ptr);
	png_destroy_write_struct(&png_ptr,(png_infopp)NULL);
	fclose(fp);
	
	return (0);
} 
