#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/profiler.h>
#include <flint/fmpz_mod_mpoly.h>
#include <flint/gr.h>

int is_perfect_power(const fmpz_t n) {
    fmpz_t temp;
    fmpz_init(temp);

    int result;
    result = fmpz_is_perfect_power(temp, n);
    fmpz_clear(temp);
    return result;
}

void get_order(fmpz *res, const fmpz *x, const fmpz *m) {
    fmpz_t i, tmp;
    fmpz_init_set_ui(i, 1);
    fmpz_init(tmp);

    while (fmpz_cmp(i, m) < 0) {
        fmpz_powm(tmp, x, i, m);
        
        if (fmpz_is_one(tmp)) {
            fmpz_set(res, i);
            fmpz_clear(i);
            fmpz_clear(tmp);
            return;
        }
        fmpz_add_ui(i, i, 1);
    }
    
    fmpz_zero(res);
    fmpz_clear(i);
    fmpz_clear(tmp);
}

void get_bound(fmpz *bound, const fmpz *r, ulong log_2_n) {
    fmpz_t phi_r;
    fmpz_init(phi_r);
    
    fmpz_euler_phi(phi_r, r);
    fmpz_sqrt(bound, phi_r);
    fmpz_mul_ui(bound, bound, log_2_n);
    fmpz_add_ui(bound, bound, 1); // because loop is excluisve
    
    fmpz_clear(phi_r);
}

int calculate_polynomials(const fmpz *n, const fmpz *r, const fmpz *bound) {
    fmpz_mod_mpoly_ctx_t ctx;
    fmpz_mod_mpoly_ctx_init(ctx, 1, ORD_LEX, n);
    
    fmpz_t n_mod_r;
    fmpz_init(n_mod_r);
    fmpz_mod(n_mod_r, n, r);

    ulong const exponent_zero[] = { 0 };
    ulong const exponent_one[] = { 1 };
    fmpz const *exponent_n_mod[] = { n_mod_r };
    fmpz const *exponent_r[] = { r };

    fmpz_mod_mpoly_t modulus, p, q, placeholder;

    // modulus: X^r - 1
    fmpz_mod_mpoly_init(modulus, ctx);
    fmpz_mod_mpoly_set_coeff_si_ui(modulus, -1, exponent_zero, ctx);
    fmpz_mod_mpoly_set_coeff_ui_fmpz(modulus, 1, (fmpz *const *)exponent_r, ctx);
    
    fmpz_t a;
    fmpz_init(a);
    for (fmpz_one(a); fmpz_cmp(a, bound) < 0; fmpz_add_ui(a, a, 1)) {
        // p = (X + a)^n
        fmpz_mod_mpoly_init(p, ctx);
        fmpz_mod_mpoly_set_coeff_fmpz_ui(p, a, exponent_zero, ctx);
        fmpz_mod_mpoly_set_coeff_ui_ui(p, 1, exponent_one, ctx);
        fmpz_mod_mpoly_pow_fmpz(p, p, n, ctx);
        
        // p mod modulus
        fmpz_mod_mpoly_init(placeholder, ctx);
        fmpz_mod_mpoly_divrem(placeholder, p, p, modulus, ctx);
        
        // q = X^n + a mod modulus
        fmpz_mod_mpoly_init(q, ctx);
        fmpz_mod_mpoly_set_coeff_fmpz_ui(q, a, exponent_zero, ctx);
        fmpz_mod_mpoly_set_coeff_ui_fmpz(q, 1, (fmpz *const *)exponent_n_mod, ctx);

        if (fmpz_mod_mpoly_cmp(p, q, ctx) != 0) {
            fmpz_mod_mpoly_clear(p, ctx);
            fmpz_mod_mpoly_clear(q, ctx);
            fmpz_mod_mpoly_clear(modulus, ctx);
            fmpz_mod_mpoly_clear(placeholder, ctx);
            fmpz_mod_mpoly_ctx_clear(ctx);
            fmpz_clear(a);

            return 0;
        }
    }

    fmpz_mod_mpoly_clear(p, ctx);
    fmpz_mod_mpoly_clear(q, ctx);
    fmpz_mod_mpoly_clear(modulus, ctx);
    fmpz_mod_mpoly_clear(placeholder, ctx);
    fmpz_mod_mpoly_ctx_clear(ctx);
    fmpz_clear(a);
    
    return 1;
}

int fixed_r_aks_test(const fmpz *r, const fmpz *n) {
    // 1. If (n = a^b for a ∈ N and b > 1), output COMPOSITE.
    if (is_perfect_power(n)) {
        return 0;
    }
    
    // 3. If 1 < gcd(a, n) < n for some a ≤ r, output COMPOSITE.
    fmpz_t gcd_result, a;
    fmpz_init(gcd_result);
    fmpz_init(a);

    for (fmpz_set_ui(a, 1); fmpz_cmp(a, r) <= 0; fmpz_add_ui(a, a, 1)) {
        fmpz_gcd(gcd_result, a, n);

        // 1 < gcd(a, n) < n
        if (fmpz_cmp_ui(gcd_result, 1) > 0 && fmpz_cmp(gcd_result, n) < 0) {            fmpz_clear(r);
            fmpz_clear(a);
            fmpz_clear(gcd_result);
            
            return 0;
        }
    }
    
    fmpz_clear(gcd_result);
    fmpz_clear(a);
    
    // 4. If n ≤ r, output PRIME.
    if (fmpz_cmp(n, r) <= 0)
    {
        return 1;   
    }
    
    // 5. For a = 1 to floor(sqrt(φ(r)) * log n) do
    //      if ((X + a)^n != X^n + a (mod X^r − 1, n)), output COMPOSITE.
    fmpz_t bound;
    fmpz_init(bound);
    ulong log_2_n = fmpz_flog_ui(n, 2); // maybe clog = ceil?
    get_bound(bound, r, log_2_n);

    int result = calculate_polynomials(n, r, bound);
    
    fmpz_clear(bound);

    // 6. Output PRIME.
    return result;
}

int is_perfect_power(const fmpz_t n) {
    fmpz_t temp;
    fmpz_init(temp);

    int result;
    result = fmpz_is_perfect_power(temp, n);
    fmpz_clear(temp);
    return result;
}

void get_order(fmpz *res, const fmpz *x, const fmpz *m) {
    fmpz_t i, tmp;
    fmpz_init_set_ui(i, 1);
    fmpz_init(tmp);

    while (fmpz_cmp(i, m) < 0) {
        fmpz_powm(tmp, x, i, m);
        
        if (fmpz_is_one(tmp)) {
            fmpz_set(res, i);
            fmpz_clear(i);
            fmpz_clear(tmp);
            return;
        }
        fmpz_add_ui(i, i, 1);
    }
    
    fmpz_zero(res);
    fmpz_clear(i);
    fmpz_clear(tmp);
}

void get_bound(fmpz *bound, const fmpz *r, ulong log_2_n) {
    fmpz_t phi_r;
    fmpz_init(phi_r);
    
    fmpz_euler_phi(phi_r, r);
    fmpz_sqrt(bound, phi_r);
    fmpz_mul_ui(bound, bound, log_2_n);
    fmpz_add_ui(bound, bound, 1); // because loop is excluisve
    
    fmpz_clear(phi_r);
}

int calculate_polynomials(const fmpz *n, const fmpz *r, const fmpz *bound) {
    fmpz_mod_mpoly_ctx_t ctx;
    fmpz_mod_mpoly_ctx_init(ctx, 1, ORD_LEX, n);
    
    fmpz_t n_mod_r;
    fmpz_init(n_mod_r);
    fmpz_mod(n_mod_r, n, r);

    ulong const exponent_zero[] = { 0 };
    ulong const exponent_one[] = { 1 };
    fmpz const *exponent_n_mod[] = { n_mod_r };
    fmpz const *exponent_r[] = { r };

    fmpz_mod_mpoly_t modulus, p, q, placeholder;

    // modulus: X^r - 1
    fmpz_mod_mpoly_init(modulus, ctx);
    fmpz_mod_mpoly_set_coeff_si_ui(modulus, -1, exponent_zero, ctx);
    fmpz_mod_mpoly_set_coeff_ui_fmpz(modulus, 1, (fmpz *const *)exponent_r, ctx);
    
    fmpz_t a;
    fmpz_init(a);
    for (fmpz_one(a); fmpz_cmp(a, bound) < 0; fmpz_add_ui(a, a, 1)) {
        // p = (X + a)^n
        fmpz_mod_mpoly_init(p, ctx);
        fmpz_mod_mpoly_set_coeff_fmpz_ui(p, a, exponent_zero, ctx);
        fmpz_mod_mpoly_set_coeff_ui_ui(p, 1, exponent_one, ctx);
        fmpz_mod_mpoly_pow_fmpz(p, p, n, ctx);
        
        // p mod modulus
        fmpz_mod_mpoly_init(placeholder, ctx);
        fmpz_mod_mpoly_divrem(placeholder, p, p, modulus, ctx);
        
        // q = X^n + a mod modulus
        fmpz_mod_mpoly_init(q, ctx);
        fmpz_mod_mpoly_set_coeff_fmpz_ui(q, a, exponent_zero, ctx);
        fmpz_mod_mpoly_set_coeff_ui_fmpz(q, 1, (fmpz *const *)exponent_n_mod, ctx);

        if (fmpz_mod_mpoly_cmp(p, q, ctx) != 0) {
            fmpz_mod_mpoly_clear(p, ctx);
            fmpz_mod_mpoly_clear(q, ctx);
            fmpz_mod_mpoly_clear(modulus, ctx);
            fmpz_mod_mpoly_clear(placeholder, ctx);
            fmpz_mod_mpoly_ctx_clear(ctx);
            fmpz_clear(a);

            return 0;
        }
    }

    fmpz_mod_mpoly_clear(p, ctx);
    fmpz_mod_mpoly_clear(q, ctx);
    fmpz_mod_mpoly_clear(modulus, ctx);
    fmpz_mod_mpoly_clear(placeholder, ctx);
    fmpz_mod_mpoly_ctx_clear(ctx);
    fmpz_clear(a);
    
    return 1;
}

int fixed_r_aks_test(const fmpz *r, const fmpz *n) {
    // 1. If (n = a^b for a ∈ N and b > 1), output COMPOSITE.
    if (is_perfect_power(n)) {
        return 0;
    }
    
    // 3. If 1 < gcd(a, n) < n for some a ≤ r, output COMPOSITE.
    fmpz_t gcd_result, a;
    fmpz_init(gcd_result);
    fmpz_init(a);

    for (fmpz_set_ui(a, 1); fmpz_cmp(a, r) <= 0; fmpz_add_ui(a, a, 1)) {
        fmpz_gcd(gcd_result, a, n);

        // 1 < gcd(a, n) < n
        if (fmpz_cmp_ui(gcd_result, 1) > 0 && fmpz_cmp(gcd_result, n) < 0) {            fmpz_clear(r);
            fmpz_clear(a);
            fmpz_clear(gcd_result);
            
            return 0;
        }
    }
    
    fmpz_clear(gcd_result);
    fmpz_clear(a);
    
    // 4. If n ≤ r, output PRIME.
    if (fmpz_cmp(n, r) <= 0)
    {
        return 1;   
    }
    
    // 5. For a = 1 to floor(sqrt(φ(r)) * log n) do
    //      if ((X + a)^n != X^n + a (mod X^r − 1, n)), output COMPOSITE.
    fmpz_t bound;
    fmpz_init(bound);
    ulong log_2_n = fmpz_flog_ui(n, 2); // maybe clog = ceil?
    get_bound(bound, r, log_2_n);

    int result = calculate_polynomials(n, r, bound);
    
    fmpz_clear(bound);

    // 6. Output PRIME.
    return result;
}

int main(int argc, char * argv[]) {
    fmpz_t r;
    fmpz_init(r);

    if (argc == 2) {
        // allow expression
        gr_ctx_t ZZ;
        gr_ctx_init_fmpz(ZZ);

        if (gr_set_str(r, argv[1], ZZ) != GR_SUCCESS) {
            flint_printf("unable to parse integer\n");
            return 1;
        }
        fmpz_t i;
        fmpz_init(i);

        flint_printf("\nVerifying the algorithm up to n = 1000 with r = %{fmpz}\n", r);
        TIMEIT_ONCE_START
        for (fmpz_one(i); fmpz_cmp_ui(i, 1000) <= 0; fmpz_add_ui(i, i, 1)) {
            
            if (fmpz_fdiv_ui(i, 100) == 0)
                flint_printf("%{fmpz}...\n", r);

            if (fixed_r_aks_test(r, i) != fmpz_is_prime(i))
                flint_printf("FAIL: wrong result for n = %{fmpz}\n", r);
        }
        flint_printf("OK\n");
        TIMEIT_ONCE_STOP
        fmpz_clear(i);

    flint_cleanup_master();
    }

    return 0;
}