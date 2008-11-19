#ifndef GRAPH_H
#define GRAPH_H

#include <algorithm>
#include <cassert>
#include <limits>

/*
 *  Compute a maximal independent set for a graph stored in CSR format
 *  using a greedy serial algorithm
 *
 *  Parameters
 *      num_rows   - number of rows in A (number of vertices)
 *      Ap[]       - CSR row pointer
 *      Aj[]       - CSR index array
 *      active     - value used for active vertices        (input)
 *       C         - value used to mark non-MIS vertices   (ouput)
 *       F         - value used to mark MIS vertices       (output)
 *      x[]        - state of each vertex
 *
 *  
 *  Returns:
 *      The number of nodes in the MIS.
 *
 *  Notes:
 *      Only the vertices with values with x[i] == active are considered 
 *      when determining the MIS.  Upon return, all active vertices will
 *      be assigned the value C or F depending on whether they are in the 
 *      MIS or not.
 *
 */
template<class I, class T>
I maximal_independent_set_serial(const I num_rows,
                                 const I Ap[], 
                                 const I Aj[], 
                                 const T active,
                                 const T  C,
                                 const T  F,
                                       T  x[])
{
    I N = 0;
    
    for(I i = 0; i < num_rows; i++){
        if(x[i] != active) continue;

        x[i] = C;
        N++;

        for(I jj = Ap[i]; jj < Ap[i+1]; jj++){
            const I j = Aj[jj];
            if(x[j] == active) {
                x[j] = F;
            }
        }

    }

    return N;
}

/*
 *  Compute a maximal independent set for a graph stored in CSR format
 *  using a variant of Luby's parallel MIS algorithm
 *
 *  Parameters
 *      num_rows   - number of rows in A (number of vertices)
 *      Ap[]       - CSR row pointer
 *      Aj[]       - CSR index array
 *      active     - value used for active vertices        (input)
 *       C         - value used to mark non-MIS vertices   (ouput)
 *       F         - value used to mark MIS vertices       (output)
 *      x[]        - state of each vertex
 *      y[]        - random values for each vertex
 *      max_iters  - maximum number of iterations 
 *                   by default max_iters=-1 and no limit 
 *                   is imposed
 *  
 *  Returns:
 *      The number of nodes in the MIS.
 *
 *  Notes:
 *      Only the vertices with values with x[i] == active are considered 
 *      when determining the MIS.  Upon return, all active vertices will
 *      be assigned the value C or F depending on whether they are in the 
 *      MIS or not.
 *  
 */
template<class I, class T, class R>
I maximal_independent_set_parallel(const I num_rows,
                                   const I Ap[], 
                                   const I Aj[],
                                   const T active,
                                   const T  C,
                                   const T  F,
                                         T  x[],
                                   const R  y[],
                                   const I  max_iters=-1)
{
    I N = 0;
    I num_iters = 0;

    bool active_nodes = true;

    while(active_nodes && (max_iters == -1 || num_iters < max_iters)){
        active_nodes = false;

        num_iters++;
        
        for(I i = 0; i < num_rows; i++){
            const R yi = y[i];

            if(x[i] != active) continue;
            
            const I row_start = Ap[i];
            const I row_end   = Ap[i+1];
    
            I jj;

            for(jj = row_start; jj < row_end; jj++){
                const I j  = Aj[jj];
                const T xj = x[j];

                if(xj == C) {
                    x[i] = F;                      //neighbor is MIS
                    break;  
                }
                
                if(xj == active){
                    const R yj = y[j];
                    if(yj > yi)
                        break;                     //neighbor is larger 
                    else if (yj == yi && j > i)
                        break;                     //tie breaker goes to neighbor
                }
            }
   
            if(jj == row_end){
                for(jj = row_start; jj < row_end; jj++){
                    const I j  = Aj[jj];
                    if(x[j] == active)
                        x[j] = F;
                }
                N++;
                x[i] = C;
            } else {
                active_nodes = true;
            }
        }
    } // end while
        
    //std::cout << std::endl << "Luby's finished in " << num_iters << " iterations " << std::endl;

    return N;
}

/*
 *  Compute a vertex coloring for a graph stored in CSR format.
 *
 *  The coloring is computed by removing maximal independent sets
 *  of vertices from the graph.
 *
 *  Returns the K, the number of colors used in the coloring.
 *  On return x[i] \in [0,1, ..., K - 1] will contain the color
 *  of the i-th vertex.
 *
 */
template<class I, class T>
T vertex_coloring_mis(const I num_rows,
                      const I Ap[], 
                      const I Aj[], 
                            T  x[])
{
    std::fill( x, x + num_rows, -1);

    I N = 0;
    T K = 0;

    while(N < num_rows){
        N += maximal_independent_set_serial(num_rows,Ap,Aj,-1-K,K,-2-K,x);
        K++;
    }

    return K;
}

template<class I, class T>
void vertex_coloring_first_fit(const I num_rows,
                               const I Ap[], 
                               const I Aj[], 
                                     T  x[],
                               const T  K)
{
    for(I i = 0; i < num_rows; i++){
        if(x[i] != K) continue;
        std::vector<bool> mask(K,false);
        for(I jj = Ap[i]; jj < Ap[i+1]; jj++){
            const I j = Aj[jj];
            if(  i == j  ) continue; //ignore diagonal
            if( x[j] < 0 ) continue; //ignore uncolored vertices
            mask[x[j]] = true;
        }
        x[i] = std::find(mask.begin(), mask.end(), false) - mask.begin();            
    }
}



/*
 * Compute a vertex coloring of a graph using the Jones-Plassmann algorithm
 *
 *  Parameters   
 *      num_rows   - number of rows in A (number of vertices)
 *      Ap[]       - CSR row pointer
 *      Aj[]       - CSR index array
 *      x[]        - color of each vertex
 *      y[]        - initial random values for each vertex
 *
 *  Notes:
 *      Arrays x and y will be overwritten
 *
 *  References:
 *      Mark T. Jones and Paul E. Plassmann
 *      A Parallel Graph Coloring Heuristic
 *      SIAM Journal on Scientific Computing 14:3 (1993) 654--669
 *      http://citeseer.ist.psu.edu/jones92parallel.html
 *      
 */
template<class I, class T, class R>
T vertex_coloring_jones_plassmann(const I num_rows,
                                  const I Ap[], 
                                  const I Aj[], 
                                        T  x[],
                                        R  y[])
{
    std::fill( x, x + num_rows, -1);

    for(I i = 0; i < num_rows; i++){
        y[i] += Ap[i+1] - Ap[i];
    }

    I N = 0;
    T K = 0; //iteration number

    while(N < num_rows){
        N += maximal_independent_set_parallel(num_rows,Ap,Aj,-1,K,-2,x,y,1);
        for(I i = 0; i < num_rows; i++){
            if(x[i] == -2)
                x[i] = -1;
        }
        vertex_coloring_first_fit(num_rows,Ap,Aj,x,K);
        K++;
    }

    return *std::max_element(x, x + num_rows);
}


/*
 * Compute a vertex coloring of a graph using the parallel 
 * Largest-Degree-First (LDF) algorithm
 *
 *  Parameters   
 *      num_rows   - number of rows in A (number of vertices)
 *      Ap[]       - CSR row pointer
 *      Aj[]       - CSR index array
 *      x[]        - color of each vertex
 *      y[]        - initial random values for each vertex
 *
 *   References:
 *     TODO 
 */
template<class I, class T, class R>
T vertex_coloring_LDF(const I num_rows,
                      const I Ap[], 
                      const I Aj[], 
                            T  x[],
                      const R  y[])
{
    std::fill( x, x + num_rows, -1);

    std::vector<R> weights(num_rows);

    I N = 0;
    T K = 0; //iteration number

    while(N < num_rows){
        // weight is # edges in induced subgraph + random value
        for(I i = 0; i < num_rows; i++){
            if(x[i] != -1) continue;
            I num_neighbors = 0;
            for(I jj = Ap[i]; jj < Ap[i+1]; jj++){
                I j = Aj[jj];
                if(x[j] == -1 && i != j)
                    num_neighbors++;
            }
            weights[i] = y[i] + num_neighbors;
        }

        N += maximal_independent_set_parallel(num_rows,Ap,Aj,-1,K,-2,x,&weights[0],1);
        for(I i = 0; i < num_rows; i++){
            if(x[i] == -2)
                x[i] = -1;
        }
        vertex_coloring_first_fit(num_rows,Ap,Aj,x,K);
        K++;
    }

    return *std::max_element(x, x + num_rows);
}

    
template<class I, class T>
void bellman_ford(const I num_rows,
                  const I Ap[], 
                  const I Aj[], 
                  const T Ax[],
                        T  x[],
                        I  y[])
{
    for(I i = 0; i < num_rows; i++){
        T xi = x[i];
        I yi = y[i];
        for(I jj = Ap[i]; jj < Ap[i+1]; jj++){
            const I j = Aj[jj];
            const T d = Ax[jj] + x[j];
            if(d < xi){
                xi = d;
                yi = y[j];
            }
        }
        x[i] = xi;
        y[i] = yi;
    }
}

// x[num_rows]     - distance to nearest seed
// y[num_rows]     - cluster membership
// z[num_centers]  - cluster centers
template<class I, class T>
void lloyd_cluster(const I num_rows,
                   const I Ap[], 
                   const I Aj[], 
                   const T Ax[],
                   const I num_seeds,
                         T  x[],
                         I  y[],
                         I  z[])
{
    for(I i = 0; i < num_rows; i++){
        x[i] = std::numeric_limits<T>::max();
        y[i] = -1;
    }
    for(I i = 0; i < num_seeds; i++){
        I seed = z[i];
        assert(seed >= 0 && seed < num_rows);
        x[seed] = 0;
        y[seed] = i;
    }

    std::vector<T> old_distances(num_rows);

    // propagate distances outward
    do{
        std::copy(x, x+num_rows, old_distances.begin());
        bellman_ford(num_rows, Ap, Aj, Ax, x, y);
    } while ( !std::equal( x, x+num_rows, old_distances.begin() ) );

    //find boundaries
    for(I i = 0; i < num_rows; i++){
        x[i] = std::numeric_limits<T>::max();
    }
    for(I i = 0; i < num_rows; i++){
        for(I jj = Ap[i]; jj < Ap[i+1]; jj++){
            I j = Aj[jj];
            if( y[i] != y[j] ){
                x[i] = 0;
                break;
            }
        }
    }

    // propagate distances inward
    do{
        std::copy(x, x+num_rows, old_distances.begin());
        bellman_ford(num_rows, Ap, Aj, Ax, x, y);
    } while ( !std::equal( x, x+num_rows, old_distances.begin() ) );


    // compute new seeds
    for(I i = 0; i < num_rows; i++){
        const I seed = y[i];

        if (seed == -1) //node belongs to no cluster
            continue;
        
        assert(seed >= 0 && seed < num_seeds);

        if( x[z[seed]] < x[i] )
            z[seed] = i;
    }
}



template<typename IndexType, typename ValueType>
void csr_propagate_max(const IndexType  num_rows,
                       const IndexType  Ap[], 
                       const IndexType  Aj[],
                       const IndexType  i_keys[],
                             IndexType  o_keys[],
                       const ValueType  i_vals[],
                             ValueType  o_vals[])
{
    for(IndexType i = 0; i < num_rows; i++){

        IndexType k_max = i_keys[i];
        ValueType v_max = i_vals[i];

        for(IndexType jj = Ap[i]; jj < Ap[i+1]; jj++){
            const IndexType j   = Aj[jj];
            const IndexType k_j = i_keys[j];
            const ValueType v_j = i_vals[j];

            if( k_j == k_max ) continue;
            if( v_j < v_max ) continue;
            if( v_j > v_max || k_j > k_max ){
                k_max = k_j;
                v_max = v_j;
            }
        }

        o_keys[i] = k_max;
        o_vals[i] = v_max;
    }
}


template<class I, class T, class R>
void maximal_independent_set_k_parallel(const I num_rows,
                                        const I Ap[], 
                                        const I Aj[],
                                        const I  k,
                                              T  x[],
                                        const R  y[],
                                        const I  max_iters=-1)
{
    std::vector<bool> active(num_rows,true);

    std::vector<I> i_keys(num_rows);
    std::vector<I> o_keys(num_rows);
    std::vector<R> i_vals(num_rows); 
    std::vector<R> o_vals(num_rows); 

    for(I i = 0; i < num_rows; i++){
        i_keys[i] = i;
        i_vals[i] = y[i];
        x[i] = 0;
    }

    for(I iter = 0; max_iters == -1 || iter < max_iters; iter++){
        for(I i = 0; i < k; i++){
            csr_propagate_max(num_rows, Ap, Aj, &(i_keys[0]), &(o_keys[0]), &(i_vals[0]), &(o_vals[0]));
            std::swap(i_keys, o_keys);
            std::swap(i_vals, o_vals);
        }

        for(I i = 0; i < num_rows; i++){
            if( i_keys[i] == i && active[i]){
                x[i] = 1; // i is a MIS-k node
            } 
            
            i_keys[i] = i;
            i_vals[i] = x[i];
        }
       
        I rank = 0;
        //while(rank < k && 2*(k - rank) > k){
        //    csr_propagate_max(num_rows, Ap, Aj, &(i_keys[0]), &(o_keys[0]), &(i_vals[0]), &(o_vals[0]));
        //    std::swap(i_keys, o_keys);
        //    std::swap(i_vals, o_vals);
        //    rank++;
        //}
        
        while(rank < k){
            csr_propagate_max(num_rows, Ap, Aj, &(i_keys[0]), &(o_keys[0]), &(i_vals[0]), &(o_vals[0]));
            std::swap(i_keys, o_keys);
            std::swap(i_vals, o_vals);
            rank++;
        }

        bool work_left = false;

        for(I i = 0; i < num_rows; i++){
            if(i_vals[i] == 1){
                active[i] =  false;
                i_vals[i] = -1;
            } else {
                i_vals[i] = y[i];
                work_left = true;
            } 
            i_keys[i] = i;
        }

        if( !work_left )
            return;
    }

}


// level must be initialized to -1!
template <class I>
void breadth_first_search(const I Ap[], 
                          const I Aj[],
                          const I seed,
                                I order[],
                                I level[])
{
    // initialize seed
    order[0]    = seed;
    level[seed] = 0;
   
    I N = 1;
    I level_begin = 0;
    I level_end   = N;

    I current_level = 1;

    while(level_begin < level_end){
        // for each node of the last level
        for(I ii = level_begin; ii < level_end; ii++){
            const I i = order[ii];

            // add all unmarked neighbors to the queue
            for(I jj = Ap[i]; jj < Ap[i+1]; jj++){
                const I j = Aj[jj];
                if(level[j] == -1){
                    order[N] = j;
                    level[j] = current_level;
                    N++;
                }
            }

        }

        level_begin = level_end;
        level_end   = N;
        current_level++;
    }

}


template <class I>
void connected_components_helper(const I Ap[], 
                                 const I Aj[],
                                 const I i,
                                 const I component,
                                       I components[])
{
    components[i] = component;

    for(I jj = Ap[i]; jj < Ap[i+1]; jj++){
        const I j = Aj[jj];
        if(components[j] == -1){
            connected_components_helper(Ap, Aj, j, component, components);
        }
    }
}

template <class I>
void connected_components(const I num_nodes,
                          const I Ap[], 
                          const I Aj[],
                                I components[])
{
    std::fill(components, components + num_nodes, -1);

    I component = 0;

    for(I i = 0; i < num_nodes; i++){
        if(components[i] == -1){
            connected_components_helper(Ap, Aj, i, component, components);
            component++;
        }
    }

}

#endif

