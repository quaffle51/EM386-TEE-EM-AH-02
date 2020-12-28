## -*- encoding: utf-8 -*-
## This file (sage_test.sagetex.sage) was *autogenerated* from sage_test.tex with sagetex.sty version 2020/08/12 v3.5.
import sagetex
_st_ = sagetex.SageTeXProcessor('sage_test', version='2020/08/12 v3.5', version_check=True)
try:
 _st_.current_tex_line = 8
 _st_.inline(0, latex(number_of_partitions(1269)))
except:
 _st_.goboom(8)
_st_.current_tex_line = 14
_st_.blockbegin()
try:
     f(x) = exp(x) * sin(2*x)
except:
 _st_.goboom(16)
_st_.blockend()
try:
 _st_.current_tex_line = 21
 _st_.inline(1, latex(f(x)))
except:
 _st_.goboom(21)
try:
 _st_.current_tex_line = 22
 _st_.inline(2, latex(diff(f, x, 2)(x)))
except:
 _st_.goboom(22)
try:
 _st_.current_tex_line = 27
 _st_.plot(0, format='notprovided', _p_=plot(f, -1, 1))
except:
 _st_.goboom(27)
_st_.endofdoc()
