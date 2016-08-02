#include "testing.h"
#include "crumple.c"
#ifdef _MPI
#include "parallel_crumple.c"
#endif

int tests_run = 0;

_test(allocation) {
    state_stack *a = new_stack("AGCU");
    assert(a);
    free_stack(a);
}_tset


_test(new_state){
    state_stack *a = new_stack("AGCU");
    new_state(a);
    state s = a->states[a->current];
    assert(s.structure[2] == -2, 
           "Structure improperly created. %d%d%d%d", 
           s.structure[0],
           s.structure[1],
           s.structure[2],
           s.structure[3]);

    assert(s.intervals[s.current_interval].end == 4, 
           "interval improperly created. %d", 
           s.intervals[s.current_interval].end);
    free_stack(a);
}_tset


_test(duplicate_state){
    state_stack *a = new_stack("AGCU");
    new_state(a);
    duplicate_state(&a->states[0], a);
    assert(a->current == 1, "No new state created!");
    state s = a->states[a->current];
    assert(s.structure[2] == -2, 
           "Structure improperly copied. %d%d%d%d", 
           s.structure[0],
           s.structure[1],
           s.structure[2],
           s.structure[3]);

    assert(s.intervals[s.current_interval].end == 4, 
           "intervals improperly copied. %d", 
           s.intervals[s.current_interval].end);
    free_stack(a);
}_tset


_test(print_state){
    state_stack *a = new_stack("AGCU");
    new_state(a);
    assert(a->states[0].intervals[0].end == 4);
    free_stack(a);
}_tset


_test(fold1){
    state_stack *a = new_stack("GCUCUAAAAGAGAG");
    new_state(a);
    OPTION.print = 0;
    while(a->current >= 0){
        fold1(a);
    }
    assert(OPTION.count == 119);
    OPTION.print = 1;
    OPTION.count = 0;

    free_stack(a);
}_tset


_test(can_pair){
    assert(can_pair('A', 'U') == 1);
    assert(can_pair('U', 'A') == 1);

    assert(can_pair('G', 'C') == 1);
    assert(can_pair('C', 'G') == 1);

    assert(can_pair('G', 'U') == 1);
    assert(can_pair('U', 'G') == 1);

    assert(can_pair('A', 'C') == 0);
    assert(can_pair('A', 'G') == 0);

    OPTION.noGU = 1;

    assert(can_pair('U', 'G') == 0);
    OPTION.noGU = 0;
    
}_tset


_test(death_triad){

    state_stack *a = new_stack("AAACCCUUU");
    OPTION.constraints = 1; 
    a->constraints = interpret_constraints(".M.......");
    new_state(a);
    make_pair(&a->states[a->current], 0, 8);
    make_pair(&a->states[a->current], 1, 7);
    make_pair(&a->states[a->current], 2, 6);
    assert(death_triad(a, 1) == TRUE, "If all three are paired and none are GU, it must be a death triad.");
    free_stack(a);


    OPTION.constraints = 0; 
    a = new_stack("AAACCCUUU");
    OPTION.constraints = 1; 
    a->constraints = interpret_constraints(".M.......");
    new_state(a);
    make_pair(&a->states[a->current], 1, 7);
    make_pair(&a->states[a->current], 2, 6);

    assert(death_triad(a, 1) == FALSE, "Both neighbors must be paired to make a death triad." );
    free_stack(a);


    OPTION.constraints = 0; 
    a = new_stack("GAACCCUUU");
    OPTION.constraints = 1; 
    a->constraints = interpret_constraints(".M.......");
    new_state(a);
    make_pair(&a->states[a->current], 0, 8);
    make_pair(&a->states[a->current], 1, 7);
    make_pair(&a->states[a->current], 2, 6);

    assert(death_triad(a, 1) == FALSE, "GU pairs should deny death triads.");
    free_stack(a);
    OPTION.constraints = 0; 
}_tset
 
_test(chemical_modification){
    state_stack *a = new_stack("AAACCCUUU");
    OPTION.constraints = 1; 
    a->constraints = interpret_constraints(".M.......");
    new_state(a);
    make_pair(&a->states[a->current], 0, 8);
    make_pair(&a->states[a->current], 1, 7);
    assert(chemical_modification_ok(a, 2, 6) == FALSE);
    free_stack(a);

    OPTION.constraints = 0; 
    a = new_stack("AAACCCUUU");
    OPTION.constraints = 1; 
    a->constraints = interpret_constraints(".M.......");
    new_state(a);
    assert(chemical_modification_ok(a, 2, 6) == TRUE);
    free_stack(a);
    OPTION.constraints = 0; 
}_tset

_test(not_pair){
    state_stack *a = new_stack("AAACCCUUU");
    OPTION.constraints = 1; 
    a->constraints = interpret_constraints("..X......");
    new_state(a);
    make_pair(&a->states[a->current], 0, 8);
    make_pair(&a->states[a->current], 1, 7);
    assert(constrained_to_not_pair(a, 2, 6) == TRUE);
    free_stack(a);

    OPTION.constraints = 0; 
    a = new_stack("AAACCCUUU");
    OPTION.constraints = 1; 
    a->constraints = interpret_constraints(".........");
    new_state(a);
    assert(constrained_to_not_pair(a, 2, 6) == FALSE);
    free_stack(a);
    OPTION.constraints = 0; 
}_tset

_test(lonely_pair){
    // This test is inadequate
    state_stack *a = new_stack("GCUCUAAAAGAGAG");
    new_state(a);
    OPTION.noLP = 1;
    OPTION.print = 0;
    while(a->current >= 0){
        fold1(a);
    }
    assert(OPTION.count == 67);
    OPTION.print = 1;
    OPTION.noLP = 0;
    OPTION.count = 0;

    free_stack(a);
}_tset

_test(helix_filter){
    // This test is inadequate.
    state_stack *a = new_stack("GCUCUAAAAGAGAG");
    new_state(a);
    OPTION.helix_length = 4;
    OPTION.helix_max_mismatch = 1;
    OPTION.helix_max_bulge = 0;
    OPTION.helix_min_count = 1;
    OPTION.print = 0;
    while(a->current >= 0){
        fold1(a);
    }
    OPTION.helix_length = 0;
    OPTION.helix_max_mismatch = 0;
    OPTION.helix_max_bulge = 0;
    OPTION.helix_min_count = 0;

    assert(OPTION.count == 10);
    OPTION.print = 1;
    OPTION.noLP = 0;
    OPTION.count = 0;

    free_stack(a);
}_tset


#ifdef _MPI

_test(split){
    OPTION.allpair = 1;
    OPTION.print = 0;
    OPTION.count = 0;
    OPTION.breadcrumb = 0;
    state_stack *a = new_stack("GCGCGCGCGCGGCGCG");
    new_state(a);
    fold1(a);
    state *b = split_stack(a);
    //    print_state(a, 0);
    //    print_stack(a);
    while(a->current >= 0){       
        //print_stack(a);
        fold1(a);
    }
    long ret1 = OPTION.count;
    //assert(OPTION.count == 161);

    free_stack(a);

    OPTION.count = 0;
                   
    a = new_stack("GCGCGCGCGCGGCGCG");
    duplicate_state(b, a);
    fold1(a);
    state *c = split_stack(a);
    //    print_state(a, 0);
    while(a->current >= 0){       
        //print_stack(a);
        fold1(a);
    }
    long ret2 = OPTION.count;

    free_stack(a);



    OPTION.count = 0;
                   
    a = new_stack("GCGCGCGCGCGGCGCG");
    duplicate_state(c, a);
    //    print_state(a, 0);
    while(a->current >= 0){
        //print_stack(a);
        fold1(a);
    }
    long ret3 = OPTION.count;

    free_stack(a);



    OPTION.count = 0;
                   
    a = new_stack("GCGCGCGCGCGGCGCG");

    new_state(a);
    while(a->current >= 0){
        fold1(a);
    }

    assert(ret1 + ret2 + ret3 == OPTION.count, 
           "%ld, %ld, %ld: %ld != %ld",
           ret1, ret2, ret3, ret1 + ret2 + ret3,  OPTION.count);

    free_stack(a);
    free_state(b);
    free_state(c);
    OPTION.print = 1;
    OPTION.count = 0;
    OPTION.allpair = 0;


}_tset


_test(pack_unpack){
    state_stack *a = new_stack("GCGCGCGCGCGGCGCG");
    new_state(a);
    state *c = split_stack(a);
    int *m = pack_state(c, a->length);
    assert(m, "Failed to pack...");
    state *d = unpack_state(m);
    assert(m, "Failed to unpack...");
    assert(c->current_interval == d->current_interval);
    assert(c->current_base == d->current_base);
    assert(c->unpair == d->unpair);
    assert(c->max_base == d->max_base);
    assert(c->intervals[c->current_interval].start == 
           d->intervals[d->current_interval].start);
    assert(c->structure[4] == 
           d->structure[4]);
} _tset

#endif

int main(int argc, char *argv[]){
    //    printf("%d, %d", );
    get_args(argc, argv);
    run_test(allocation);
    run_test(new_state);
    run_test(duplicate_state);
    run_test(print_state);
    run_test(fold1);
    run_test(can_pair);
    run_test(death_triad);
    run_test(chemical_modification);
    run_test(not_pair);
    run_test(lonely_pair);
    run_test(helix_filter);
#ifdef _MPI
   run_test(split);
   run_test(pack_unpack);
#endif
    printf("PASSED\n");
    return 0;
}
