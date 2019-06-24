#   ZONOTOPE LIBRARY ver 1.2
#
#   Copyright (c) 2017, Marco Cococcioni.
#   All rights reserved.
#
#   Redistribution and use in source and binary forms, with or without
#   modification, are permitted provided that the following conditions are met:
#
#   1. Redistributions of source code must retain the above copyright notice, this
#      list of conditions and the following disclaimer.
#   2. Redistributions in binary form must reproduce the above copyright notice,
#      this list of conditions and the following disclaimer in the documentation
#      and/or other materials provided with the distribution.
#
#   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#   ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#   WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#   DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
#   ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
#   (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
#   LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
#   ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

zonotope <- function(star, computeAdditionalStatistics) {
  
  if (missing(computeAdditionalStatistics))
    computeAdditionalStatistics <- T
  
  dim <- ncol(star)   # number of dimensions
  nn <- nrow(star)    # number of generators
  
  if (computeAdditionalStatistics){
      printf("-----------------------------------------------------------------------");
      printf("ZONOTOPE LIBRARY VER 1.2                                               ");
      printf("-----------------------------------------------------------------------");	
      printf("INPUT: SET OF GENERATORS                                               ");
      printf("N. of dimensions (including the output): %d                            ", dim );
      printf("N. of generators: %d                                                   ", nrow(star) );
      printf("-----------------------------------------------------------------------");
      printf("... the computation of the volume has started (it can take a while) ...");
  }
  
  # Start the clock!
  ptm <- proc.time()
  
  if ( nn < dim ){
    print('Zonotope::Error: the number of generators must be greater than the dimension of the space');
    return;
  }
  
  # compute_vertices <- FALSE
  # tagi <- matrix(0,nn,1)
  # # compute tagi
  
  face <- matrix(0, nn, 1)
  
  volume <- 0
  
  counters <- matrix(0, dim, 1)
  
  for (i in 1:dim - 1) {
    counters[i] <- i
  }
  counters[dim] <- nn + 1
  
  a <- matrix(0, dim, dim)
  
  
  issue_linear_dependency_warning <- FALSE
  keep_going <- TRUE
  
  while (keep_going) {
    a[1:(dim-1), (1:dim)] <- star[counters[1:(dim-1)], 1:dim]
    a[dim, 1:dim] <- star[counters[1], 1:dim]
    
    # Compute the normal of the vectors in a
    # A normal can be obtained by computing the null space of a (i.e., its kernel).
    
    normal <- Null(t(a)) # normal <- kernel(a)
    
    if (ncol(normal) > 1){
        if ( issue_linear_dependency_warning){
            print("WARNING! Some of the generators are linearly dependent. The volume might be inaccurate")
            issue_linear_dependency_warning <- F
        }  
        normal <- t(t(normal[, 1])) # take first column only
    }
    
    normal <- normal / norm(normal, type="F")
    
    # I need to find orientation of the normal
    a[dim, 1:dim] <- normal[1:dim]
    s <- det(a)

    if (s < 0)
      normal <- -normal
    
    aux <- sign(star %*% normal)
    
    starting <- 1
    
    for (i in 1:dim) {
      
      ending <- counters[i]
      
      if (ending-1 >= starting){
          K <- starting:ending-1
          face[K] <- aux[K]
      }
      if (ending < nn) {
        face[ending] <- 0
      }
      starting <- counters[i] + 1
    }
    
    base <- colSums( star * face[,rep(1,dim)] )
    
    pyramid <- crossprod(base, normal * abs(s)) # scalar product between base (a row vector)
                                                # and normal*abs(s), a column vector
    volume <- volume + pyramid

    # increase counters
    aux <- increase( counters, dim - 2, dim - 1, nrow(star) )
     
    keep_going <- aux$keep_going
    
    counters   <- aux$counters # rimettere sopra questo assegnamento!!!
  
  } # while (1)
  
  volume <- volume / dim

  # Stop the clock
  etime <- proc.time() - ptm

  printf("-----------------------------------------------------------------------");
  printf("OUTPUT: STATISTICS                                                     ");
  printf("S1: Total volume:  %g                                                  ", volume);
  
  if (computeAdditionalStatistics){
  
      # computing diagonal and total_length_squared
      # diagonal <- matrix(0,1,dim)
      diagonal <- colSums(star)
      
      total_length_squared <- 0
      for (i in 1:nn)
           total_length_squared = total_length_squared + sum(star[i,]^2)
      cat("S2: Diagonal: ", sprintf('%.4f ',diagonal), "\n");
   
      norm_of_the_diagonal <- sqrt(sum(diagonal^2))
      printf("S3: Diagonal norm: %g", norm_of_the_diagonal)
  
      printf("S4: Sum of squared norms: %g", total_length_squared)

      # Compute the Gini index (it depends on volume and diagonal)
      gini <- volume / prod(diagonal);
      printf("S5: Gini index: %g", gini);
     
      base <- 0
      for ( i in 1: (dim - 1) )
          base = base + diagonal[i]*diagonal[i]
      tang_diag_input <- diagonal[dim] / sqrt(base)
     
      printf("S6: Tangent of angle btw. diagonal and the input plane: %g", tang_diag_input)
     
      cos1 = diagonal[dim] / norm_of_the_diagonal
      printf("S7: Cosine against output: %g", cos1 )
  
      cos2 = diagonal[1]/sqrt(sum((diagonal[1:(dim-1)])^2))
        
      printf("S8: Cosine of projection of diagonal on input plane with x axis: %g", cos2)
     
      printf("S9: Volume against the cube of the norm of the diagonal: %g" , volume/( norm_of_the_diagonal*norm_of_the_diagonal*norm_of_the_diagonal ) );
     
      #        aux = (total_length_squared/3.0).^1.5;
      #        mistery_number = binomial(size(star,1),3) * aux;
      #        cat("B6: Mystery number: %g\n", mistery_number);       
      #        cat("B7: Tangent against input axes: %g\n", diagonal(2)/diagonal(1) );               
      #        cat("B8: Tangent against input axes: %g\n", xx);        
      #        cat("B9: Solid angle: %g\n", solidAngle(star, edges));        
      #        cat("B10: Normalized vectors volume: %g\n", xxx);                
      #        cat("B11: Volume against diagonal cubed of boundary vectors: %g \n", xxxx);        
  
  }
  
  printf("-----------------------------------------------------------------------");
  cat("Elapsed time (MIN): \n");
  printf("%g", etime[3]/60);    
  printf("-----------------------------------------------------------------------");
  
  return(volume)
  
} # End of the ZONOTOPE function


# -----------------------------------------------------------
# UTILITY FUNCTIONS
# -----------------------------------------------------------
increase <- function(counters, pos, dim, nn){

  if ( pos < 0 ){
    keep_going <- FALSE
    return(list(keep_going=keep_going,counters=counters))
  }


  if ( counters[pos+1] == ( nn - (dim - pos) + 1 ) ){
    aux <- increase(counters, pos-1, dim, nn)
    keep_going <- aux$keep_going
    counters   <- aux$counters
    return(list(keep_going=keep_going,counters=counters))
  }

  counters[pos+1] <- counters[pos+1] + 1

  if ( (pos+1) <= (dim-1) ){
    for ( i in (pos+1):(dim-1) )
      counters[i+1] <- counters[(pos+1)]+(i-pos)
  }
  keep_going <- TRUE
  return (list(keep_going=keep_going, counters=counters))
}


Null<-function(M)
{
  tmp <- qr(M)
  set <- if(tmp$rank == 0L) seq_len(ncol(M)) else  -seq_len(tmp$rank)
  qr.Q(tmp, complete = TRUE)[, set, drop = FALSE]
}

printf<- function(s, ...){
  cat(paste0(sprintf(s, ...)), '\n')
}

#------------------------------------------------
# UNUSED FUNCTIONS
#------------------------------------------------

# myseq<-function(from,to){
#  if (from <= to){
#    return(seq(from,to,1))
#  }else
#    return(NULL)
# }


# repmat <- function(X,m,nn){
#   ##R equivalent of repmat (matlab)
#   mx = dim(X)[1]
#   nx = dim(X)[2]
#   return(matrix(t(matrix(X,mx,nx*nn)),mx*m,nx*nn,byrow=T))
# }



# # [ project, release 2.0 ]
# project <- function(face, star, dim){
#   aux = repmat(face,1,dim)
#   vND = colSums( star * aux )
# #  vND = 1.200549  4.019081  6.231622  4.794844 -2.640008     in demo1
#   return(vND)
# }

# [ project, release 1.0 ]
# project <- function(face, star){  
#   dim = ncol(star)
#   v3D <- matrix(0,1,dim)
#   nn <- nrows(star)
#   for ( i in 1:nn ){
#       v3D <- v3D + Vec3_optimesdouble( star[i,], face[i,] )
#   }
#   return(V3D)
# }
# 
# Vec3_optimesdouble <- function (a3D, doublescalar){
#   v3D = a3D*doublescalar
#   return(V3D)
# 
# }

# normalize <- function(v){
#   normalized_v = v / norm(v, type="F")
#   return(normalized_v)
# }


# % function newvertices = find4vertices(face,star,i,j)
# % dim = 3;
# % newvertices = zeros(4,dim);
# % 
# % face(i) = -1;
# % face(j) = -1;  %-1, -1
# % newvertices(1,:) = project(face, star);
# % 
# % face(i) = 1;   %1, -1
# % newvertices(2,:) = project(face, star);
# % 
# % face(j) = 1;   %1, 1
# % newvertices(3,:) = project(face, star);
# % 
# % face(i) = -1;   %-1 1
# % newvertices(4,:) = project(face, star);
# % end % find4vertices
#     
#     
