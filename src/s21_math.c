#include "s21_math.h"

#include <stdio.h>

int s21_abs(int x) {
  int res = 0;
  if (x >= 0) {
    res = x;
  } else {
    res = -1 * x;
  }
  return res;
}

long double s21_log(double x) {
  double y_n = 0, y_n1 = 0;
  if (x > 0) {
    while (1) {
      y_n1 = 2 * (x - s21_exp(y_n)) / (x + s21_exp(y_n));
      if (s21_fabs(y_n1) > S21_EPS) {
        y_n += y_n1;
      } else {
        break;
      }
    }
    if (x == S21_POS_INF) {
      y_n = S21_POS_INF;
    }
  } else if (x == 0) {
    y_n = S21_NEG_INF;
  } else {
    y_n = S21_NAN;
  }
  return (long double)y_n;
}

long double s21_exp(double x) {
  long double out = 0;
  long double addentum = 0;
  unsigned long long n = 1;
  double fact = 1;
  double power = 1;
  if ((x != S21_POS_INF) && (x != S21_NEG_INF)) {
    while (1) {
      addentum = (long double)power / fact;
      if (s21_fabs((double)addentum) > S21_EPS) {
        out += addentum;
      } else {
        break;
      }
      power = x * power;
      fact *= n;
      n++;
      if (out == S21_POS_INF || out == S21_NEG_INF) {
        break;
      }
    }
  } else if (x == S21_POS_INF) {
    out = S21_POS_INF;
  } else if (x == S21_NEG_INF) {
    out = 0;
  }
  return out;
}

long double s21_pow(double base, double exp) {
  int sign = 1;
  long double out = s21_fabs(base);
  int check_fractional = (int)(exp * 10) % 10;
  int check_even = (int)(exp * 10) / 10;
  if ((base < 0) && (check_fractional != 0)) {
    out = S21_NAN;
    if (base == S21_NEG_INF) {
      out = 0;
      if (exp > 0) {
        out = S21_POS_INF;
      }
    } else if (exp == S21_NEG_INF) {
      out = 0.0;
    } else if (exp == S21_POS_INF) {
      out = S21_POS_INF;
    }
  } else {
    if (base == 0) {
      out = 0.0;
      if (exp == 0) {
        out = 1.0;
      } else if (exp < 0) {
        out = S21_POS_INF;
      }
    } else {
      if ((base < 0) && (check_even % 2 == 0) && (check_fractional == 0)) {
        sign = -1;
        base = sign * base;
        out = s21_exp(exp * s21_log(base));
      } else {
        sign = 1;
        out = sign * s21_exp(exp * s21_log(base));
      }
      if ((base < 0) && (check_even % 2 == 1) && (check_fractional == 0)) {
        sign = -1;
        base = sign * base;
        out = sign * s21_exp(exp * s21_log(base));
      } else {
        sign = 1;
        out = sign * s21_exp(exp * s21_log(base));
      }
    }
  }
  return out;
}

long double s21_acos(double x) {
  long double out = 0;
  double addentum = x;
  double n = 1;
  double step = 0;
  if (s21_fabs(x) < 1) {
    while (1) {
      step = x * x * (2 * n - 1) * (2 * n - 1) / (4 * n + 2) / n;
      if (s21_fabs(addentum) > S21_EPS * 10) {
        out += (long double)addentum;
      } else {
        break;
      }
      addentum *= step;
      n++;
    }
    out = S21_PI_2 - out;
  } else if (x == -1) {
    out = S21_PI;
  } else if (x == 1) {
    out = 0;
  } else {
    out = 0.0 / 0.0;
  }
  return out;
}

long double s21_fabs(double x) {
  long double res = x;
  if (res < 0) {
    res *= -1;
  }
  return res;
}

long double s21_sin(double x) {
  long double res = 0;
  long double addentum = 0;
  long double n = 0;
  int sign = 1;
  if ((x != S21_POS_INF) && (x != S21_NEG_INF)) {
    x = s21_fmod(x, 2 * S21_PI);
    while (1) {
      addentum = sign * s21_pow(x, (2 * n + 1)) / s21_factorial(2 * n + 1);
      if (s21_fabs(addentum) > S21_EPS) {
        res += addentum;
      } else {
        break;
      }
      sign *= -1;
      n++;
    }
  } else {
    res = S21_NAN;
  }
  return res;
}

long double s21_cos(double x) {
  long double res = 0;
  long double addentum = 0;
  long long n = 0;
  int sign = 1;
  if ((x != S21_POS_INF) && (x != S21_NEG_INF)) {
    x = s21_fmod(x, 2 * S21_PI);

    while (1) {
      addentum = sign * s21_pow(x, (2 * n)) / s21_factorial(2 * n);
      if (s21_fabs(addentum) > S21_EPS * 100) {
        res += addentum;
      } else {
        break;
      }
      sign *= -1;
      n += 1;
    }
  } else {
    res = S21_NAN;
  }
  return res;
}

long double s21_tan(double x) {
  long double out = 0;
  if ((x != S21_POS_INF) && (x != S21_NEG_INF)) {
    if (x == S21_PI_2) {
      out = S21_POS_INF;
    }
    if (x == -S21_PI_2) {
      out = S21_NEG_INF;
    }
    if (s21_fabs(x) != S21_PI_2) {
      x = s21_fmod(x, S21_PI_2);
      out = s21_sin(x) / s21_cos(x);
    }
  } else {
    out = S21_NAN;
  }
  return out;
}

long double s21_ceil(double x) {
  long double out = x;
  long long int_num = (long long)x;
  long long int_num_plus = int_num + 1;
  long double double_num = (long double)int_num;
  long double double_num_plus = (long double)int_num_plus;
  if ((x != S21_POS_INF) && (x != S21_NEG_INF)) {
    if (x < 0) {
      out = double_num;
    } else {
      if (out != double_num) {
        out = double_num_plus;
      } else {
        out = double_num;
      }
    }
  } else if (x == S21_POS_INF) {
    out = S21_POS_INF;
  } else if (x == S21_NEG_INF) {
    out = S21_NEG_INF;
  }
  return out;
}

long double s21_floor(double x) {
  long double out = x;
  long long int_num = (long long)x;
  long long int_num_minus = int_num - 1;
  long double double_num = (long double)int_num;
  long double double_num_minus = (long double)int_num_minus;
  if ((x != S21_POS_INF) && (x != S21_NEG_INF)) {
    if (x > 0) {
      out = double_num;
    } else {
      if (out != double_num) {
        out = double_num_minus;
      } else {
        out = double_num;
      }
    }
  } else if (x == S21_POS_INF) {
    out = S21_POS_INF;
  } else if (x == S21_NEG_INF) {
    out = S21_NEG_INF;
  }
  return out;
}

long double s21_factorial(long double x) {
  long double out = 1;
  if (x <= 1) {
    out = 1;
  } else {
    for (int i = 1; i <= x; i++) {
      out *= i;
    }
  }
  return out;
}

long double s21_fmod(double x, double y) {
  long double res = 0;
  if ((y != S21_POS_INF) && (y != S21_NEG_INF) && (x != S21_POS_INF) &&
      (x != S21_NEG_INF)) {
    res = (x / y - (int)(x / y)) * y;
  } else {
    if (x == 0) {
      res = 0;
    } else {
      res = x;
      if ((x == S21_POS_INF) || (x == S21_NEG_INF)) {
        res = S21_NAN;
      }
    }
  }

  return res;
}

long double s21_sqrt(double x) {
  long double left = 0;
  long double right = s21_right_border(x);
  long double middle = 0;
  long double res = 0;
  if (x > 0) {
    middle = left + (right - left) / 2;
    while (s21_fabs(middle * middle - x) > S21_EPS * 10000) {
      if (middle * middle < x) {
        left = middle;
      } else {
        right = middle;
      }
      middle = left + (right - left) / 2;
    }
    res = middle;
  } else if (x == 0) {
    res = 0;
  } else if (x < 0) {
    res = S21_NAN;
  }
  if (x == S21_POS_INF) {
    res = S21_POS_INF;
  }
  return res;
}

long double s21_asin(double x) {
  long double res = 0;
  long double addentum = 0;
  long double n = 0;
  if (x < 1 && x > -1 && x != 0) {
    while (1) {
      addentum =
          (s21_pow(x, (2 * n + 1)) * s21_factorial(2 * n)) /
          (s21_pow(4, n) * s21_factorial(n) * s21_factorial(n) * (2 * n + 1));
      if (s21_fabs(addentum) > S21_EPS * 1E+6) {
        res += addentum;
      } else {
        break;
      }
      n += 1;
    }
  } else if (x == 1) {
    res = S21_PI / 2;
  } else if (x == -1) {
    res = -S21_PI / 2;
  } else if (x == 0) {
    res = 0.0;
  } else {
    res = S21_NAN;
  }
  return res;
}

long double s21_right_border(long double x) {
  long double res = 0;
  if (x > 1) {
    res = x;
  } else {
    res = 1;
  }
  return res;
}

long double s21_atan(double x) {
  long double res = 0;
  if ((x != S21_POS_INF) && (x != S21_NEG_INF)) {
    res = s21_asin(x / s21_sqrt(1 + x * x));
  } else {
    if (x == S21_POS_INF) {
      res = S21_PI_2;
    } else if (x == S21_NEG_INF) {
      res = -S21_PI_2;
    }
  }
  return res;
}
