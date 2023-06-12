While installing an error appeared in this file with:

In file included from /usr/lib64/gcc/x86_64-suse-linux/13/include/immintrin.h:57,
                 from /home/xxx/anaconda3/lib/python3.9/site-packages/numpy/distutils/checks/cpu_avx512_knl.c:14:
In function ‘_mm512_mask_prefetch_i64scatter_pd’,
    inlined from ‘main’ at /home/jean/anaconda3/lib/python3.9/site-packages/numpy/distutils/checks/cpu_avx512_knl.c:23:5:
/usr/lib64/gcc/x86_64-suse-linux/13/include/avx512pfintrin.h:180:3: error: ‘base’ may be used uninitialized [-Werror=maybe-uninitialized]
  180 |   __builtin_ia32_scatterpfqpd (__mask, (__v8di) __index, __addr, __scale,
      |   ^~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  181 |                                __hint);
      |                                ~~~~~~~
<built-in>: In function ‘main’:
<built-in>: note: by argument 3 of type ‘const void *’ to ‘__builtin_ia32_scatterpfqpd’ declared here
/home/jean/anaconda3/lib/python3.9/site-packages/numpy/distutils/checks/cpu_avx512_knl.c:18:9: note: ‘base’ declared here
   18 |     int base[128];
      |         ^~~~
cc1: all warnings being treated as errors

In order to "fix" this I used:

#include <immintrin.h>

int main(int argc, char **argv)
{
    int base[128] = {0};  // Initialize the base array with zeros
    __m512d ad = _mm512_loadu_pd((const __m512d*)argv[argc-1]);
    /* ER */
    __m512i a = _mm512_castpd_si512(_mm512_exp2a23_pd(ad));
    /* PF */
    _mm512_mask_prefetch_i64scatter_pd(base, _mm512_cmpeq_epi64_mask(a, a), a, 1, _MM_HINT_T1);
    return base[0];
}

which is different from the original file attached. Why is this error happening?
