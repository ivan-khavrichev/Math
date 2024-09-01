all: s21_math.a

WFLAGS = -Wall -Werror -Wextra
TEST_FLAGS = -lcheck
CC = gcc

s21_math.a:
	gcc $(WFLAGS) -c "s21_math.c" "s21_math.h"
	ar rc s21_math.a s21_math.o
	ranlib s21_math.a

test:
	gcc $(WFLAGS) $(TEST_FLAGS) s21_math.c unit_tests.c -o unit_tests

gcov_report:
	$(CC) --coverage unit_tests.c s21_math.c $(WFLAGS) $(TEST_FLAGS) -o test
	./test
	lcov -t "test" -o test.info -c -d .
	genhtml -o report test.info
	open report/index.html

clean:
	rm -rf *.o *.g* *.a *.out s21_math test.info test unit_tests report

rebuild: clean all