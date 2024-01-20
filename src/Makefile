CC=gcc
CFLAGS=-Wall -Werror -Wextra -std=c11
VALGRIND=valgrind --leak-check=full --leak-resolution=med -s
OS:= $(shell uname -s)
LINFLAGS=-lcheck -pthread -lrt -lm -lsubunit

SRCFILE = s21_matrix.c
OBJECTS = s21_matrix.o
TESTS =  s21_test.c

all: clean test

s21_matrix.a:
	$(CC) -c $(SRCFILE)
	@ar rcs s21_matrix.a $(OBJECTS)
	ranlib s21_matrix.a

test: s21_matrix.a
ifeq ($(OS), Darwin)
	gcc $(TESTS) s21_matrix.a  -o s21_test $(shell pkg-config --cflags --libs check) -lcheck -lpthread #-fsanitize=address
	./s21_test
	rm -rf *.o *.a
else
	gcc $(TESTS) s21_matrix.a -o s21_test $(LINFLAGS)
	./s21_test
	rm -rf *.o *.a
endif

gcov_report: clean
ifeq ($(OS), Darwin)
	gcc --coverage $(TESTS) -o s21_test $(SRCFILE) $(shell pkg-config --cflags --libs check) -lcheck
	./s21_test
	gcov -f s21_matrix.gcda
	lcov -t "s21_test" -o s21_test.info -c -d .   
	genhtml -o report s21_test.info
	open ./report/index.html
else
	# sudo apt-get install xdg-utils
	gcc --coverage $(TESTS) -o s21_test $(SRCFILE) $(LINFLAGS)
	./s21_test
	gcov -f s21_test-s21_matrix.gcda
	lcov -t "s21_test" -o s21_test.info -c -d .  
	genhtml -o report s21_test.info
	xdg-open ./report/index.html
endif

check: clean cppcheck style leaks

cppcheck: test
	cppcheck --enable=all --suppress=missingIncludeSystem *.c *.h

style: test
	cp ../materials/linters/.clang-format .clang-format 
	clang-format -style=google -i *.c *.h
	clang-format -style=google -n *.c *.h
	rm .clang-format 

leaks: test
ifeq ($(OS), Darwin)
	CK_FORK=no leaks --atExit -- ./s21_test
	rm -rf s21_matrix.o
else
	# sudo apt install -y valgrind
	$(VALGRIND) ./s21_test
	rm -rf s21_matrix.o
endif

clean:
	-rm -rf *.o && rm -rf *.gcno && rm -rf *.exec
	-rm -rf *.a && rm -rf *.gcda
	-rm -rf *.info && rm -rf *.gcov
	-rm -rf ./gcov_report
	-rm -rf ./run_tests ./s21_test 
	-rm -rf ./report/
