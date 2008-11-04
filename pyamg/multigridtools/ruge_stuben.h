#ifndef RUGE_STUBEN_H
#define RUGE_STUBEN_H

#include <iostream>
#include <vector>
#include <iterator>
#include <cassert>


/*
 *
 * Return a strength of connection matrix using the classical 
 * strength of connection measure by Ruge and Stuben.
 *
 * Specifically, an off-diagonal entry A[i.j] is a strong 
 * connection if:
 *  
 *      -A[i,j] >= theta * max( -A[i,k] )   where k != i
 * 
 * Otherwise, the connection is weak.
 *  
 */          

template<class I, class T, class F>
void classical_strength_of_connection(const I n_row,
                                      const F theta,
                                      const I Ap[], const I Aj[], const T Ax[],
                                            I Sp[],       I Sj[],       T Sx[])
{
    I nnz = 0;
    Sp[0] = 0;

    for(I i = 0; i < n_row; i++){
        F max_offdiagonal = std::numeric_limits<F>::min();

        const I row_start = Ap[i];
        const I row_end   = Ap[i+1];

        for(I jj = row_start; jj < row_end; jj++){
            if(Aj[jj] != i){
                max_offdiagonal = std::max(max_offdiagonal,mynorm(Ax[jj]));
            }
        }

        F threshold = theta*max_offdiagonal;
        for(I jj = row_start; jj < row_end; jj++){
            F norm_jj = mynorm(Ax[jj]);
            if(norm_jj >= threshold){
                if(Aj[jj] != i){
                    Sj[nnz] = Aj[jj];
                    Sx[nnz] = Ax[jj];
                    nnz++;
                }
            }
        }

        Sp[i+1] = nnz;
    }
}


#define F_NODE 0
#define C_NODE 1
#define U_NODE 2

template<class I>
void rs_cf_splitting(const I n_nodes,
                     const I Sp[], const I Sj[], 
                     const I Tp[], const I Tj[], 
                           I splitting[])
{
    std::vector<I> lambda(n_nodes,0);

    //compute lambdas
    for(I i = 0; i < n_nodes; i++){
        lambda[i] = Tp[i+1] - Tp[i];
    }

    //for each value of lambda, create an interval of nodes with that value
    // ptr - is the first index of the interval
    // count - is the number of indices in that interval
    // index to node - the node located at a given index
    // node to index - the index of a given node
    std::vector<I> interval_ptr(n_nodes+1,0);
    std::vector<I> interval_count(n_nodes+1,0);
    std::vector<I> index_to_node(n_nodes);
    std::vector<I> node_to_index(n_nodes);

    for(I i = 0; i < n_nodes; i++){
        interval_count[lambda[i]]++;
    }
    for(I i = 0, cumsum = 0; i < n_nodes; i++){
        interval_ptr[i] = cumsum;
        cumsum += interval_count[i];
        interval_count[i] = 0;
    }
    for(I i = 0; i < n_nodes; i++){
        I lambda_i = lambda[i];
        I index    = interval_ptr[lambda_i] + interval_count[lambda_i];
        index_to_node[index] = i;
        node_to_index[i]     = index;
        interval_count[lambda_i]++;
    }


    std::fill(splitting, splitting + n_nodes, U_NODE);

    //Now add elements to C and F, in decending order of lambda
    for(I top_index = n_nodes - 1; top_index != -1; top_index--){
        I i        = index_to_node[top_index];
        I lambda_i = lambda[i];

        //if (n_nodes == 4)
        //    std::cout << "selecting node #" << i << " with lambda " << lambda[i] << std::endl;

        //remove i from its interval
        interval_count[lambda_i]--;

        if(splitting[i] == F_NODE)
        {
            continue;
        } 
        else 
        {
            assert(splitting[i] == U_NODE); 

            splitting[i] = C_NODE;

            //For each j in S^T_i /\ U
            for(I jj = Tp[i]; jj < Tp[i+1]; jj++){
                I j = Tj[jj];

                if(splitting[j] == U_NODE){
                    splitting[j] = F_NODE;

                    //For each k in S_j /\ U
                    for(I kk = Sp[j]; kk < Sp[j+1]; kk++){
                        I k = Sj[kk];

                        if(splitting[k] == U_NODE){	      
                            //move k to the end of its current interval
                            if(lambda[k] >= n_nodes - 1) continue;

                            // TODO make this robust
                            //if(lambda[k] >= n_nodes -1)
                            //    std::cout << std::endl << "lambda[" << k << "]=" << lambda[k] << " n_nodes=" << n_nodes << std::endl;
                            //assert(lambda[k] < n_nodes - 1);//this would cause problems!
                            
                            I lambda_k = lambda[k];
                            I old_pos  = node_to_index[k];
                            I new_pos  = interval_ptr[lambda_k] + interval_count[lambda_k] - 1;

                            node_to_index[index_to_node[old_pos]] = new_pos;
                            node_to_index[index_to_node[new_pos]] = old_pos;
                            std::swap(index_to_node[old_pos], index_to_node[new_pos]);

                            //update intervals
                            interval_count[lambda_k]   -= 1;
                            interval_count[lambda_k+1] += 1; //invalid write!
                            interval_ptr[lambda_k+1]    = new_pos;

                            //increment lambda_k
                            lambda[k]++;
                        }
                    }
                }
            }

            //For each j in S_i /\ U
            for(I jj = Sp[i]; jj < Sp[i+1]; jj++){
                I j = Sj[jj];
                if(splitting[j] == U_NODE){            //decrement lambda for node j
                    if(lambda[j] == 0) continue;

                    //assert(lambda[j] > 0);//this would cause problems!

                    //move j to the beginning of its current interval
                    I lambda_j = lambda[j];
                    I old_pos  = node_to_index[j];
                    I new_pos  = interval_ptr[lambda_j]; 

                    node_to_index[index_to_node[old_pos]] = new_pos;
                    node_to_index[index_to_node[new_pos]] = old_pos;
                    std::swap(index_to_node[old_pos],index_to_node[new_pos]);

                    //update intervals
                    interval_count[lambda_j]   -= 1;
                    interval_count[lambda_j-1] += 1;
                    interval_ptr[lambda_j]     += 1;
                    interval_ptr[lambda_j-1]    = interval_ptr[lambda_j] - interval_count[lambda_j-1];

                    //decrement lambda_j
                    lambda[j]--;
                }
            }
        }
    }
}

/*
 *   Produce the Ruge-Stuben prolongator using "Direct Interpolation"
 *
 *
 *   The first pass uses the strength of connection matrix 'S' 
 *   and C/F splitting to compute the row pointer for the prolongator.
 *
 *   The second pass fills in the nonzero entries of the prolongator
 *
 *   Reference:
 *      Page 479 of "Multigrid"
 *
 */      
template<class I>
void rs_direct_interpolation_pass1(const I n_nodes,
                                   const I Sp[], const I Sj[],
                                   const I splitting[],
                                         I Bp[])
{
    I nnz = 0;
    Bp[0] = 0;
    for(I i = 0; i < n_nodes; i++){
        if( splitting[i] == C_NODE ){
            nnz++;
        } else {
            for(I jj = Sp[i]; jj < Sp[i+1]; jj++){
                if (splitting[Sj[jj]] == C_NODE)
                    nnz++;
            }
        }
        Bp[i+1] = nnz;
    }
}


template<class I, class T>
void rs_direct_interpolation_pass2(const I n_nodes,
                                   const I Ap[], const I Aj[], const T Ax[],
                                   const I Sp[], const I Sj[], const T Sx[],
                                   const I splitting[],
                                   const I Bp[],       I Bj[],       T Bx[])
{

    for(I i = 0; i < n_nodes; i++){
        if(splitting[i] == C_NODE){
            Bj[Bp[i]] = i;
            Bx[Bp[i]] = 1;
        } else {
            T sum_strong_pos = 0, sum_strong_neg = 0;
            for(I jj = Sp[i]; jj < Sp[i+1]; jj++){
                if (splitting[Sj[jj]] == C_NODE){
                    if (Sx[jj] < 0)
                        sum_strong_neg += Sx[jj];
                    else
                        sum_strong_pos += Sx[jj];
                }
            }

            T sum_all_pos = 0, sum_all_neg = 0;
            T diag = 0;
            for(I jj = Ap[i]; jj < Ap[i+1]; jj++){
                if (Aj[jj] == i){
                    diag += Ax[jj];
                } else {
                    if (Ax[jj] < 0)
                        sum_all_neg += Ax[jj];
                    else
                        sum_all_pos += Ax[jj];
                }
            }

            T alpha = sum_all_neg / sum_strong_neg;
            T beta  = sum_all_pos / sum_strong_pos;

            if (sum_strong_pos == 0){
                diag += sum_all_pos;
                beta = 0;
            }

            T neg_coeff = -alpha/diag;
            T pos_coeff = -beta/diag;

            I nnz = Bp[i];
            for(I jj = Sp[i]; jj < Sp[i+1]; jj++){
                if (splitting[Sj[jj]] == C_NODE){
                    Bj[nnz] = Sj[jj];
                    if (Sx[jj] < 0)
                        Bx[nnz] = neg_coeff * Sx[jj];
                    else
                        Bx[nnz] = pos_coeff * Sx[jj];
                    nnz++;
                }
            }
        }
    }


    std::vector<I> map(n_nodes);
    for(I i = 0, sum = 0; i < n_nodes; i++){
        map[i]  = sum;
        sum    += splitting[i];
    }
    for(I i = 0; i < Bp[n_nodes]; i++){
        Bj[i] = map[Bj[i]];
    }
}




template<class I, class T>
void rs_standard_interpolation(const I n_nodes,
                               const I Ap[], const I Aj[], const T Ax[],
                               const I Sp[], const I Sj[], const T Sx[],
                               const I Tp[], const I Tj[], const T Tx[],
                                     I Bp[],       I Bj[],       T Bx[])
{

    // TODO
}


//#define NodeType char
// // The following function closely approximates the
// // method described in the 1987 Ruge-Stuben paper
//
//template<class I, class T>
//void rs_interpolation(const I n_nodes,
//        const I Ap[], const I Aj[], const T Ax[],
//        const I Sp[], const I Sj[], const T Sx[],
//        const I Tp[], const I Tj[], const T Tx[],
//        std::vector<I> * Bp, std::vector<I> * Bj, std::vector<T> * Bx){
//
//    std::vector<I> lambda(n_nodes,0);
//
//    //compute lambdas
//    for(I i = 0; i < n_nodes; i++){
//        lambda[i] = Tp[i+1] - Tp[i];
//    }
//
//
//    //for each value of lambda, create an interval of nodes with that value
//    // ptr - is the first index of the interval
//    // count - is the number of indices in that interval
//    // index to node - the node located at a given index
//    // node to index - the index of a given node
//    std::vector<I> interval_ptr(n_nodes,0);
//    std::vector<I> interval_count(n_nodes,0);
//    std::vector<I> index_to_node(n_nodes);
//    std::vector<I> node_to_index(n_nodes);
//
//    for(I i = 0; i < n_nodes; i++){
//        interval_count[lambda[i]]++;
//    }
//    for(I i = 0, cumsum = 0; i < n_nodes; i++){
//        interval_ptr[i] = cumsum;
//        cumsum += interval_count[i];
//        interval_count[i] = 0;
//    }
//    for(I i = 0; i < n_nodes; i++){
//        I lambda_i = lambda[i];
//        I index    = interval_ptr[lambda_i]+interval_count[lambda_i];
//        index_to_node[index] = i;
//        node_to_index[i]     = index;
//        interval_count[lambda_i]++;
//    }
//
//
//
//
//
//    std::vector<NodeType> NodeSets(n_nodes,U_NODE);
//
//    //Now add elements to C and F, in decending order of lambda
//    for(I top_index = n_nodes - 1; top_index > -1; top_index--){
//        I i        = index_to_node[top_index];
//        I lambda_i = lambda[i];
//#ifdef DEBUG
//        {
//#ifdef DEBUG_PRINT
//            std::cout << "top_index " << top_index << std::endl;
//            std::cout << "i         " << i << std::endl;
//            std::cout << "lambda_i  " << lambda_i << std::endl;
//
//            for(I i = 0; i < n_nodes; i++){
//                std::cout << i << "=";
//                if(NodeSets[i] == U_NODE)
//                    std::cout << "U";
//                else if(NodeSets[i] == F_NODE)
//                    std::cout << "F";
//                else
//                    std::cout << "C";
//                std::cout << " ";
//            }
//            std::cout << std::endl;
//
//            std::cout << "node_to_index" << std::endl;
//            for(I i = 0; i < n_nodes; i++){
//                std::cout << i << "->" << node_to_index[i] << "  ";
//            }
//            std::cout << std::endl;
//            std::cout << "index_to_node" << std::endl;
//            for(I i = 0; i < n_nodes; i++){
//                std::cout << i << "->" << index_to_node[i] << "  ";
//            }
//            std::cout << std::endl;
//
//            std::cout << "interval_count ";
//            for(I i = 0; i < n_nodes; i++){
//                std::cout << interval_count[i] << " ";
//            }
//            std::cout << std::endl;
//#endif
//
//            //make sure arrays are correct
//            for(I n = 0; n < n_nodes; n++){
//                assert(index_to_node[node_to_index[n]] == n);
//            }
//
//            //make sure intervals are reasonable
//            I sum_intervals = 0;
//            for(I n = 0; n < n_nodes; n++){
//                assert(interval_count[n] >= 0);
//                if(interval_count[n] > 0){
//                    assert(interval_ptr[n] == sum_intervals);
//                }
//                sum_intervals += interval_count[n];
//            }
//            assert(sum_intervals == top_index+1);
//
//
//            if(interval_count[lambda_i] <= 0){
//                std::cout << "top_index " << top_index << std::endl;
//                std::cout << "lambda_i " << lambda_i << std::endl;
//                std::cout << "interval_count[lambda_i] " << interval_count[lambda_i] << std::endl;
//                std::cout << "top_index " << top_index << std::endl;
//                std::cout << "i         " << i << std::endl;
//                std::cout << "lambda_i  " << lambda_i << std::endl;
//            }
//
//
//            for(I n = 0; n <= top_index; n++){
//                assert(NodeSets[index_to_node[n]] != C_NODE);
//            }
//        }
//        assert(node_to_index[i] == top_index);
//        assert(interval_ptr[lambda_i] + interval_count[lambda_i] - 1 == top_index);
//        //max interval should have at least one element
//        assert(interval_count[lambda_i] > 0);    
//#endif
//
//
//        //remove i from its interval
//        interval_count[lambda_i]--;
//
//
//        if(NodeSets[i] == F_NODE){
//            continue;
//        } else {
//            assert(NodeSets[i] == U_NODE);
//
//            NodeSets[i] = C_NODE;
//
//            //For each j in S^T_i /\ U
//            for(I jj = Tp[i]; jj < Tp[i+1]; jj++){
//                I j = Tj[jj];
//
//                if(NodeSets[j] == U_NODE){
//                    NodeSets[j] = F_NODE;
//
//                    //For each k in S_j /\ U
//                    for(I kk = Sp[j]; kk < Sp[j+1]; kk++){
//                        I k = Sj[kk];
//
//                        if(NodeSets[k] == U_NODE){	      
//                            //move k to the end of its current interval
//                            assert(lambda[j] < n_nodes - 1);//this would cause problems!
//
//                            I lambda_k = lambda[k];
//                            I old_pos  = node_to_index[k];
//                            I new_pos  = interval_ptr[lambda_k] + interval_count[lambda_k] - 1;
//
//                            node_to_index[index_to_node[old_pos]] = new_pos;
//                            node_to_index[index_to_node[new_pos]] = old_pos;
//                            std::swap(index_to_node[old_pos],index_to_node[new_pos]);
//
//                            //update intervals
//                            interval_count[lambda_k]   -= 1;
//                            interval_count[lambda_k+1] += 1;
//                            interval_ptr[lambda_k+1]    = new_pos;
//
//                            //increment lambda_k
//                            lambda[k]++;
//
//#ifdef DEBUG
//                            assert(interval_count[lambda_k]   >= 0);
//                            assert(interval_count[lambda_k+1] >  0);
//                            assert(interval_ptr[lambda[k]] <= node_to_index[k]);
//                            assert(node_to_index[k] < interval_ptr[lambda[k]] + interval_count[lambda[k]]);
//#endif
//                        }
//                    }
//                }
//            }
//
//            //For each j in S_i /\ U
//            for(I jj = Sp[i]; jj < Sp[i+1]; jj++){
//                I j = Sj[jj];
//                if(NodeSets[j] == U_NODE){            //decrement lambda for node j
//                    assert(lambda[j] > 0);//this would cause problems!
//
//                    //move j to the beginning of its current interval
//                    I lambda_j = lambda[j];
//                    I old_pos  = node_to_index[j];
//                    I new_pos  = interval_ptr[lambda_j]; 
//
//                    node_to_index[index_to_node[old_pos]] = new_pos;
//                    node_to_index[index_to_node[new_pos]] = old_pos;
//                    std::swap(index_to_node[old_pos],index_to_node[new_pos]);
//
//                    //update intervals
//                    interval_count[lambda_j]   -= 1;
//                    interval_count[lambda_j-1] += 1;
//                    interval_ptr[lambda_j]     += 1;
//                    interval_ptr[lambda_j-1]    = interval_ptr[lambda_j] - interval_count[lambda_j-1];
//
//                    //decrement lambda_j
//                    lambda[j]--;
//
//#ifdef DEBUG
//                    assert(interval_count[lambda_j]   >= 0);
//                    assert(interval_count[lambda_j-1] >  0);
//                    assert(interval_ptr[lambda[j]] <= node_to_index[j]);
//                    assert(node_to_index[j] < interval_ptr[lambda[j]] + interval_count[lambda[j]]);
//#endif
//                }
//            }
//        }
//    }
//
//
//
//
//#ifdef DEBUG
//    //make sure each f-node has at least one strong c-node neighbor
//    for(I i = 0; i < n_nodes; i++){
//        if(NodeSets[i] == F_NODE){
//            I row_start = Sp[i];
//            I row_end   = Sp[i+1];
//            bool has_c_neighbor = false;
//            for(I jj = row_start; jj < row_end; jj++){
//                if(NodeSets[Sj[jj]] == C_NODE){
//                    has_c_neighbor = true;
//                    break;
//                }
//            }
//            assert(has_c_neighbor);
//        }   
//    }
//#endif
//
//    //Now construct interpolation operator
//    std::vector<T> d_k(n_nodes,0);
//    std::vector<bool> C_i(n_nodes,0);
//    Bp->push_back(0);
//    for(I i = 0; i < n_nodes; i++){
//        if(NodeSets[i] == C_NODE){
//            //interpolate directly
//            Bj->push_back(i);
//            Bx->push_back(1);      
//            Bp->push_back(Bj->size());
//        } else {
//            //F_NODE
//
//            //Step 4
//            T d_i = 0; //denominator for this row
//            for(I jj = Ap[i]; jj < Ap[i+1]; jj++){ d_i += Ax[jj]; }
//            for(I jj = Sp[i]; jj < Sp[i+1]; jj++){ d_i -= Sx[jj]; }
//
//            //Create C_i, initialize d_k
//            for(I jj = Sp[i]; jj < Sp[i+1]; jj++){ 
//                I j = Sj[jj];
//                if(NodeSets[j] == C_NODE){
//                    C_i[j] = true;
//                    d_k[j] = Sx[jj];
//                }
//            }
//
//            bool Sj_intersects_Ci = true; //in the case that i has no F-neighbors
//            for(I jj = Sp[i]; jj < Sp[i+1]; jj++){ //for j in D^s_i
//                I    j = Sj[jj];
//                T   a_ij = Sx[jj];
//                T   a_jl = 0;
//
//                if(NodeSets[j] != F_NODE){continue;}
//
//                //Step 5
//                Sj_intersects_Ci = false;
//
//                //compute sum a_jl
//                for(I ll = Sp[j]; ll < Sp[j+1]; ll++){
//                    if(C_i[Sj[ll]]){
//                        Sj_intersects_Ci = true;
//                        a_jl += Sx[ll];
//                    }	    
//                }
//
//                if(!Sj_intersects_Ci){ break; }
//
//                for(I kk = Sp[j]; kk < Sp[j+1]; kk++){
//                    I   k = Sj[kk];
//                    T  a_jk = Sx[kk];
//                    if(C_i[k]){
//                        d_k[k] += a_ij*a_jk / a_jl;
//                    }	    
//                }
//            }
//
//            //Step 6
//            if(Sj_intersects_Ci){
//                for(I jj = Sp[i]; jj < Sp[i+1]; jj++){ 
//                    I j = Sj[jj];
//                    if(NodeSets[j] == C_NODE){
//                        Bj->push_back(j);
//                        Bx->push_back(-d_k[j]/d_i);      
//                    }
//                }	
//                Bp->push_back(Bj->size());
//            } else { //make i a C_NODE
//                NodeSets[i] = C_NODE;
//                Bj->push_back(i);
//                Bx->push_back(1);      
//                Bp->push_back(Bj->size());
//            }
//
//
//            //Clear C_i,d_k
//            for(I jj = Sp[i]; jj < Sp[i+1]; jj++){ 
//                I j = Sj[jj];
//                C_i[j] = false;
//                d_k[j] = 0;
//            }
//
//        }
//
//    }
//
//    //for each c-node, determine its index in the coarser lvl
//    std::vector<I> cnode_index(n_nodes,-1);
//    I n_cnodes = 0;
//    for(I i = 0; i < n_nodes; i++){
//        if(NodeSets[i] == C_NODE){
//            cnode_index[i] = n_cnodes++;
//        }
//    }
//    //map old C indices to coarse indices
//    for(typename std::vector<I>::iterator iter = Bj->begin(); iter != Bj->end(); iter++){
//        *iter = cnode_index[*iter];
//    }
//}

#endif
