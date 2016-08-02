#include <stdlib.h>
#include <stdio.h>

int len;

void usage(char *me)
{
    fprintf(stderr,
            "usage: %s infile outfile\neither or both can be - for stdin. if both are - then stdin is treated as though the files were concatenated, so the sequence should appear on the first line and the constraints are on the second line and the structures to test are on the following lines\n",
            me);
}

void skip_sequence(FILE * fp)
{
    char c;
    for (c = getc(fp); c != EOF && c != '\n'; c = getc(fp));
}

char *get_structure(FILE * fp)
{
    // TODO: check malloc/realloc return values
    char c;
    char *out;
    int incr = 2000, alloc;

    out = malloc((sizeof *out) * (alloc = incr));

    for (len = 0, c = getc(fp); c != EOF && c != '\n'; c = getc(fp)) {
        out[len++] = c;
        if (len >= alloc) {
            out = realloc(out, (sizeof *out) * (alloc += incr));
        }
    }

    out = realloc(out, (sizeof *out) * (len + 1));
    out[len] = '\0';

    return out;
}

void make_pair_table(char *C, int **pos, int **neg)
{
    int *p, *n, *ps, *ns;
    int j, ph, nh;

    p = malloc((sizeof *p) * (len + 1));
    n = malloc((sizeof *n) * (len + 1));
    ps = malloc((sizeof *ps) * (len + 1));
    ns = malloc((sizeof *ns) * (len + 1));

    ph = nh = 0;

    for (j = 1; j <= len; ++j) {
        p[j] = n[j] = 0;
        if (C[j - 1] == '(')
            ps[ph++] = j;
        if (C[j - 1] == '[')
            ns[nh++] = j;
        if (C[j - 1] == ')') {
            if (ph <= 0) {
                fprintf(stderr, "unbalanced ( ) at %d in %s\n", j, C);
                exit(1);
            }
            p[ps[--ph]] = j;
        }
        if (C[j - 1] == ']') {
            if (nh <= 0) {
                fprintf(stderr, "unbalanced [ ] at %d in %s\n", j, C);
                exit(1);
            }
            n[ns[--nh]] = j;
        }
    }

    free(ps);
    free(ns);

    if (pos)
        *pos = p;
    else
        free(p);
    if (neg)
        *neg = n;
    else
        free(n);
}

void compare(char *c, int *s, int *p, int *n)
{
    int e = 0, i, j;
    for (i = 1; i <= len; ++i) {
        j = s[i];
        if (j && n[i] == j) {
            printf("%s should not pair (%d,%d)\n", c, i, n[i]);
            e = 1;
        }
        if (p[i] && p[i] != j) {
            printf("%s should pair (%d,%d)\n", c, i, p[i]);
            e = 1;
        }
    }
    if (!e)
        printf("%s ok\n", c);
}

int main(int arfc, char **arfv)
{
    /*
       arfv[1] is the input file with sequence on the first line and constraint on the second line
       takes constraints of the ( ) and [ ] variety
       arfv[2] is the list of structures output by subopt to check
     */

    if (arfc != 3) {
        usage(arfv[0]);
        exit(1);
    }

    FILE *fin, *fout;
    int i, *pos, *neg;
    char *C;

    if (arfv[1][0] == '-' && arfv[1][1] == '\0')
        fin = stdin;
    else if (!(fin = fopen(arfv[1], "rb"))) {
        fprintf(stderr, "unable to open %s, bailing\n", arfv[1]);
        exit(1);
    }

    if (arfv[2][0] == '-' && arfv[2][1] == '\0')
        fout = stdin;
    else if (!(fout = fopen(arfv[2], "rb"))) {
        fprintf(stderr, "unable to open %s, bailing\n", arfv[2]);
        if (stdin != fin)
            fclose(fin);
        exit(1);
    }

    skip_sequence(fin);
    C = get_structure(fin);
    make_pair_table(C, &pos, &neg);

    printf("constraint string:\n%s\n", C);
    free(C);

    printf("forced pairs:");
    for (i = 1; i <= len; ++i) {
        if (pos[i])
            printf(" (%d,%d)", i, pos[i]);
    }
    printf("\n");

    printf("forbidden pairs:");
    for (i = 1; i <= len; ++i) {
        if (neg[i])
            printf(" (%d,%d)", i, neg[i]);
    }
    printf("\n");

    while (!feof(fout)) {
        char *s;
        int *ps;
        s = get_structure(fout);
        make_pair_table(s, &ps, NULL);
        if (s[0]) {
            compare(s, ps, pos, neg);
        }
        free(s);
        free(ps);
    }

    free(pos);
    free(neg);
    if (stdin != fin)
        fclose(fin);
    if (stdin != fin)
        fclose(fout);

    return 0;
}
