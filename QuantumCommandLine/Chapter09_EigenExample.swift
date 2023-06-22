/**
 
 This file contains an example for calculating eigenvalues and vectors of a  operator.
  
 There is not much documentation on how to make Swift call the LAPACK libraries that are contained in the accelerate framework (so have hardware acceleration).
 
 Here we use SGEEV that computes for an N-by-N real non symmetric matrix A the
  eigenvalues and, optionally, the left and/or right eigenvectors.
 
 We take advantage of the fact that a complex matrix can be written as a real matrix of double the size.
 
 This has only been tested on real and hermitian matrices and only just enough to reproduce a previously published figure. If this file is used as a template, heavy testing is recommended.
 
 Due to the pressure of completing the work in time for a deadline  there are numerous violations of good coding practice in this file. The most important violation here is not to fix issues as they have been seen.
 
 This file contains
    - A driver function: computeExampleSQUIDEigenvalues
    - A struct for a SQUID ring
    - A basic test of LAPAC agains a textbook example: TestEigenSolver
    - and extension to MatrixOperator that provides a wrapper function for SGEEV

 */

import Foundation
import Accelerate
import Quantum

// Reproduce figure 3 of https://arxiv.org/pdf/cond-mat/0411115.pdf
// https://doi.org/10.1103/PhysRevB.64.184517

public func Chapter10_computeExampleSQUIDEigenvalues(){
    let squidSpace = StateSpace(dimension: 20, label: "Squid")
    let squid = SQUID(space: squidSpace)
    
    let startFlux = 0.0
    let endFlux = 1.0
    let nFluxSteps = 100
    let dFlux = (endFlux - startFlux)/Double(nFluxSteps)
    
    let pathToFile = FileManager.default.homeDirectoryForCurrentUser.path + "/data/SQUID/"
    let writeFilename = pathToFile + "Eigenvalues.dat"
    _ = FileManager.default.createFile(atPath: writeFilename, contents: nil, attributes: nil)
    let fileHandle = FileHandle(forWritingAtPath: writeFilename)
    
    let headingText = "flux  ground first    second  third\n"
    fileHandle?.write(headingText.data(using: .utf8)!)
    for fluxStep in 0 ... nFluxSteps {
        let flux = startFlux + Double(fluxStep) * dFlux
        let H = squid.getSquidHamiltonian(external_flux_in_flux_quantum: flux)

        let (u,_) = H.eigensolve()
        var v = [Double]()
        for eigenvalue in u {
            v.append(eigenvalue.real/squid.hwo)
        }
        v.sort(by: <)
        let outputText = "\(flux)\t\(v[0])\t\(v[1])\t\(v[2])\t\(v[3])\n"
        fileHandle?.write(outputText.data(using: .utf8)!)
    }
}
/**
 Latex code to produce plot with pgfplot
 \begin{tikzpicture}
 \begin{axis}[   xmin=0,
                 xmax=1,
                 width=11cm,
                 xlabel = External flux $\Phi_x/\Phi_0$,
                 ylabel = Energy in $\hbar \omega_0$]
 \addplot[color=black]table[x=flux,y=ground]{Eigenvalues.dat};
 \addplot[color=black]table[x=flux,y=first]{Eigenvalues.dat};
 \addplot[color=black]table[x=flux,y=second]{Eigenvalues.dat};
 \end{axis}
 \end{tikzpicture}
 */
public struct SQUID {
    public let hwo: Double
    public let hnu: Double
    public let space: StateSpace
    private let H_0: MatrixOperator
    private var GammLn: [Double]
    private let f: Double
    private let dimension: Int
    private var NSM: MatrixOperator
    private var NCM: MatrixOperator
    // see https://arxiv.org/pdf/cond-mat/0411115.pdf for explanation of units and reference plot (Figure 3)
    public init(hwo_inPhi_oSqrOverLambda: Double = 0.043, hnu_inPhi_oSqrOverLambda: Double = 0.07, space: StateSpace) {
        self.hwo = hwo_inPhi_oSqrOverLambda
        self.hnu = hnu_inPhi_oSqrOverLambda
        self.space = space
        dimension = space.dimension

        // Store values so as to not recompute log gamma every time.
        // Old trick to deal with large factorial - can be improved upon
        GammLn = Array.init(repeating: 0.0, count: dimension)
        for i in 0 ..<  dimension {
            GammLn[i] = lgamma(Double(i+1))
        }
        // saves a small amount of computation time
        f = Double.pi * sqrt(2.0 * hwo)

        // The harmonic oscillator part of the Hamiltonian.
        H_0 =  hwo * ( space.numberOperator + 0.5 * space.identityOperator )

        NSM = space.nullOperator
        NCM = space.nullOperator
        nsm()
        ncm()
    }
    public func getSquidHamiltonian(external_flux_in_flux_quantum PhiX: Double) -> MatrixOperator {
        let bias =  2.0 * Double.pi * PhiX
        return H_0 + hnu * ( NSM * sin(bias) - NCM * cos(bias) )
    }
    /**
     * Makes a matrix of external flux independent sine and cosine terms
     *
     * Details of mathematics can be found in:
     *      Limits to the observation of coherent oscillations in a SQUID ring, Physica B 215 (1995)367 376
     *      https://doi.org/10.1016/0921-4526(95)00417-0
     *
     * Code modified from my Java/C code modified from FORTRAN code by Joe Diggins.
     *
     * Comments contain original Java code. Some original code kept as comments for traceability with java implementation.
     * 
     * Also note that as standard matrix elements these should really be extensions to in StateSpace just like the number operator.
     * This code could be substantially improved including removing the use of GammLn
     *
     */
    private mutating func nsm() {
        
        // double sum;
        // int    n,nmax;
        // double test,bd,bb;
        // int    bi,bc,row,column;

        for row in 0 ..< dimension {            // }(row=0;row<dimension;row++) {
            for column in 0 ..< dimension {     // (column=0;column<dimension;column++) {
                var sum = 0.0
                let nmax = ( (row < column) ? row : column )
                for n in 0 ... nmax {           // }(n=0;n<=nmax;++n) {
                    let bi = row + column - 2 * n
                    if bi.isMultiple(of: 2) == false { // if((bi>>1<<1) != bi) {
                        //var test = 2.0;
                        let bc = (bi-1)/2
                        var bd = f.power(UInt(bi)) * 2.0  // Math.pow(f,bi) * test;
                       // bd*=Math.exp(0.5*GammLn[row]+0.5*GammLn[column]
                       //                  -GammLn[n]-GammLn[row-n]-GammLn[column-n]);
                        let arg1 = 0.5 * ( GammLn[row] + GammLn[column] )
                        let arg2 = GammLn[n] + GammLn[row-n] + GammLn[column-n]
                        bd *= exp( arg1 - arg2 )
                        let bb = ( (bc >> 1) << 1) == bc ? 1.0 : -1.0  // -1.0^bc
                        sum += bd * bb / 2.0
                    }
                }
                NSM[row,column] = Complex(real: sum * exp(-f * f / 2.0) )
            }
        }
    }
    private mutating func ncm() {
        //    double sum;
        //    int    n,nmax;
        //    double test,bd,bb;
        //    int    bi,bc,row,column;
        
        for row in 0 ..< dimension {            // for(row=0;row<dimension;row++) {
            for column in 0 ..< dimension {     //for(column=0;column<dimension;column++) {
                var sum = 0.0
                let nmax = ( (row<column) ? row : column )
                for n in 0 ... nmax {                     // for (n=0;n<=nmax;++n) {
                    let bi = row + column - 2 * n         // Not very DRY code!
                    if bi.isMultiple(of: 2) || bi == 0 {   // if(((bi>>1<<1) == bi) || bi==0) {
                        // test = 2.0;
                        let bc = bi / 2
                        var bd = f.power(UInt(bi)) * 2.0  //  =Math.pow(f,bi)*test;
                        
                        // bd*=Math.exp(0.5*GammLn[row]+0.5*GammLn[column]
                        //  -GammLn[n]-GammLn[row-n]-GammLn[column-n]);
                        let arg1 = 0.5 * ( GammLn[row] + GammLn[column] )
                        let arg2 = GammLn[n] + GammLn[row-n] + GammLn[column-n]
                        bd *= exp( arg1 - arg2 )
                        let bb = ( (bc >> 1) << 1 ) == bc ? 1.0 : -1.0  // -1^bc
                        sum += bd * bb / 2.0
                    }
                }
                NCM[row,column] = Complex(real: sum * exp(-f * f / 2.0))
            }
        }
    }
}

// https://web.cs.ucdavis.edu/~bai/publications/baidemmeletal06.pdf page 75-13
// note in test third eigenvector is - the printed version which agrees up to a total phase
public func TestEigenSolver() {
    let space = StateSpace(dimension: 4, label: "test")
    
    let elements = [Complex(real: 4), Complex(real: -5),Complex(real: 0),Complex(real: 3),
                    Complex(real: 0), Complex(real: 4),Complex(real: -3),Complex(real: -5),
                    Complex(real: 5), Complex(real: -3),Complex(real: 4),Complex(real: 0),
                    Complex(real: 3), Complex(real: 0),Complex(real: 5),Complex(real: 4)]
    let A = MatrixOperator(elements: elements, in: space)
    let (u,v) = A.eigensolve()
    
    print(u)
    
    for vector in v {
        print(vector)
    }
}

// Code kept separate from main implementation to make clear this work was not done in the one-week of coding time frame. Should be well tested then migrated to the core library.
extension MatrixOperator {
    // see examples at
    // https://developer.apple.com/documentation/accelerate/solving_systems_of_linear_equations_with_lapack
    // A nice description of the parameters (some comments taken verbatim from this document)
    // https://web.cs.ucdavis.edu/~bai/publications/baidemmeletal06.pdf
    // Warning - only uses float rather than double precision
    public func eigensolve() -> (values: [ComplexReal], vectors: [StateVector] ) {
        let d = space.dimension

        var output_values: [ComplexReal]
        var output_vectors: [StateVector]
        output_values = Array(repeating: Complex(0.0), count: d)
        var vectorWksp = Array(repeating: Complex(0.0), count: d)
        output_vectors = []
        // is it row major or column major - needs a test
        
        
        
        // Inputs
        var JOBVL = Int8("V".utf8.first!) // calculate left eigenvector
        var JOBVR = Int8("N".utf8.first!) // do not calculate right vector
        var n = __CLPK_integer(d * 2) // order of the matrix, * 2 as complex input
        var A = Array(repeating: Float(0), count: d * d * 4 ) // the matrix of the eigen problem
        
        // store complex matrix in a larger real matrix
        let dimA = 2 * d
        for i in 0 ..< d {
            for j in 0 ..< d {
                A[  i      + dimA * j       ] = Float(self[i,j].real)
                A[ (i+d)   + dimA * (j+d)   ] = Float(self[i,j].real) // lower right = upper left
                A[  i      + dimA * (j + d) ] = Float(self[i,j].imag) // upper right
                A[ (i + d) + dimA * j       ] = Float(self[i,j].imag) // lower left = - upper right for complex arithmatic
            }
        }
        var LDA = n  // the leading dimension of the array A as square = n
        var LDVL = n // leading dimensions of the arrays VL
        var LDVR = n // leading dimensions of the arrays VL
        var WORK = Array(repeating: Float(0), count: 8 * d ) // workspace
        var LWORK = __CLPK_integer(d * 8) // dimension of workspace
        // Output
        var WR = Array(repeating: Float(0), count: d * 2 ) // real parts of the computed eigenvalues
        var WI = Array(repeating: Float(0), count: d * 2 ) // imaginary parts of the computed eigenvalues
        // Note: Complex conjugate pairs of eigenvalues appear consecutively with the eigenvalue having the positive imaginary part first.
        var VR = Array(repeating: Float(0), count: 8 * d * 8 * d ) // left eigenvectors
        var VL = Array(repeating: Float(0), count: 8 * d * 8 * d ) // right eigenvectors
        var INFO = Int32(0) // status
        
        
        var _ = sgeev_(&JOBVL, &JOBVR, &n, &A, &LDA, &WR, &WI, &VL, &LDVL, &VR, &LDVR, &WORK, &LWORK, &INFO)
        if INFO != __CLPK_integer(0) {
            //INFO: = 0 if successful exit.
            //If INFO = −j, the jth argument had an illegal value.
            //If INFO = j, the QR algorithm failed to compute all the eigenvalues, and no eigen- vectors have been computed; elements j + 1 : N of WR and WI contain eigenvalues, which have converged.

            errorStream.write("Problem in sgeev_. INFO = \(INFO)")
        }

        var i = 0
        while i < d {
            output_values[i] = Complex( real: Double(WR[i]), imag: Double(WI[i]) )
            if WI[i] == Float(0) {
                for j in 0 ..< d {
                    vectorWksp[j] = Complex(real: Double( VL[j + i * dimA] ) )
                }
                output_vectors.append(StateVector(elements: vectorWksp, in: space))
            } else { // output stored in conjugate pairs
                output_values[i+1] = Complex( real: Double(WR[i]), imag: -Double(WI[i]) )
                for j in 0 ..< d {
                    vectorWksp[j] = Complex(real: Double( VL[j + i * dimA] ),
                                            imag: Double( VL[j + (i + 1) * dimA]))
                }
                output_vectors.append(StateVector(elements: vectorWksp, in: space))
                for j in 0 ..< d {
                    vectorWksp[j] = Complex(real: Double( VL[j + i * dimA] ),
                                            imag: -Double( VL[j + (i + 1) * dimA]))
                }
                output_vectors.append(StateVector(elements: vectorWksp, in: space))
                i = i + 1
            }
            i = i + 1
        }
        
        return (output_values, output_vectors)
    }
}

//
//  Created by M J Everitt on 09/09/2022.
//
