main_iplib: main_iplib.c ip_lib.o bmp.o ip_lib.h bmp.h
	gcc main_iplib.c ip_lib.o bmp.o -o main_iplib -Wall --ansi --pedantic -lm -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra
	
ip_lib.o: ip_lib.c ip_lib.h
	gcc ip_lib.c -o ip_lib.o -c -Wall --ansi --pedantic -lm -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra
	
bmp.o: bmp.c bmp.h
	gcc bmp.c -o bmp.o -c -Wall
	
clean:
	@rm -f main_iplib ip_lib.o bmp.o
