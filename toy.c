#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define TK_K 3
#define TK_N 4
#define TK_Q 97

//polynomial multiplication in Z97[X]/(X^4+1)
static void toy_polmul_naive(short *dst, const short *a, const short *b, int add_to_dst)
{
    dst[0] = (a[0] * b[0] - a[3] * b[1] - a[2] * b[2] - a[1] * b[3]) % TK_Q;
    dst[1] = (a[1] * b[0] + a[0] * b[1] - a[3] * b[2] - a[2] * b[3]) % TK_Q;
    dst[2] = (a[2] * b[0] + a[1] * b[1] + a[0] * b[2] - a[3] * b[3]) % TK_Q;
    dst[3] = (a[3] * b[0] + a[2] * b[1] + a[1] * b[2] + a[0] * b[3]) % TK_Q;

    if (!add_to_dst)
    {
        for (int i = 0; i < TK_N; i++)
        {
            if (dst[i] < 0)
                dst[i] += TK_Q;
        }
    }
}


void toy_gen(short *A, short *t, short *s)
{
    // Fill K*K-matrix A with uniformly random numbers mod q
    for (int i = 0; i < TK_K * TK_K; i++)
    {
        A[i] = rand() % TK_Q;
    }

    // Fill K-vectors s & e with small random numbers mod q
    for (int i = 0; i < TK_K; i++)
    {
        s[i] = rand() % TK_Q;
    }

    // Fill K-vector t with zeros
    for (int i = 0; i < TK_K; i++)
    {
        t[i] = 0;
    }

    // Matrix-vector multiplication in Zq[X]/(X^n+1)
    for (int i = 0; i < TK_K; i++)
    {
        for (int j = 0; j < TK_K; j++)
        {
            short tmp[TK_N];
            toy_polmul_naive(tmp, &A[i * TK_K * TK_N + j * TK_N], &s[j * TK_N], 1);
            for (int k = 0; k < TK_N; k++)
            {
                t[i * TK_N + k] += tmp[k];
                t[i * TK_N + k] %= TK_Q;
            }
        }
    }
}

void toy_enc(const short *A, const short *t, int plain, short *u, short *v)
{
    short r[TK_K * TK_N];
    short e1[TK_K * TK_N];
    short e2[TK_N];

    // Fill K-vectors r & e1 with small random numbers mod q
    for (int i = 0; i < TK_K * TK_N; i++)
    {
        r[i] = rand() % TK_Q;
        e1[i] = rand() % TK_Q;
    }

    // Fill scalar (one-polynomial) e2 with small random numbers mod q
    for (int i = 0; i < TK_N; i++)
    {
        e2[i] = rand() % TK_Q;
    }

    // Matrix-vector multiplication in Zq[X]/(X^n+1)
    for (int i = 0; i < TK_K; i++)
    {
        u[i] = 0;
        for (int j = 0; j < TK_K; j++)
        {
            short tmp[TK_N];
            toy_polmul_naive(tmp, &A[j * TK_K * TK_N + i * TK_N], &r[j * TK_N], 1);
            for (int k = 0; k < TK_N; k++)
            {
                u[i * TK_N + k] += tmp[k];
                u[i * TK_N + k] %= TK_Q;
            }
        }
    }

    // Calculate ciphertext v
    for (int i = 0; i < TK_N; i++)
    {
        v[i] = e2[i];
        for (int j = 0; j < TK_K; j++)
        {
            v[i] += t[j * TK_N + i] * r[j * TK_N];
            v[i] %= TK_Q;
        }
        v[i] += (plain >> i) & 1 ? TK_Q / 2 : 0;
        v[i] %= TK_Q;
    }
}

int toy_dec(const short *s, const short *u, const short *v)
{
    short p[TK_N];
    int plain = 0;

    // Calculate p = v - s dot u
    for (int i = 0; i < TK_N; i++)
    {
        p[i] = v[i];
        for (int j = 0; j < TK_K; j++)
        {
            p[i] -= s[j] * u[j * TK_N + i];
            p[i] %= TK_Q;
        }
    }

    // Calculate plaintext
    for (int i = 0; i < TK_N; i++)
    {
        int val = p[i];
        if (val > TK_Q / 2)
            val -= TK_Q;
        int bit = abs(val) > TK_Q / 4;
        plain |= bit << i;
    }

    return plain;
}

int main()
{
    // Seed the random number generator
    srand(time(NULL));

    short A[TK_K * TK_K * TK_N];
    short t[TK_K * TK_N];
    short s[TK_K];
    short u[TK_K * TK_N];
    short v[TK_N];

    // Generate keys
    toy_gen(A, t, s);

    // Encrypt
    int plain = 5; // Example plaintext should be from 0 to 15
    toy_enc(A, t, plain, u, v);

    // Decrypt
    int decrypted = toy_dec(s, u, v);
    printf("Plaintext: %d\n", plain);
    printf("Decrypted: %d\n", decrypted);


    return 0;
    
}

