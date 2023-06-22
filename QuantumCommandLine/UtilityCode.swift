/*
 This file contains various utility functions that were written in support
 of each of the programmes used to generate figures for the book.
 
 Strictly speaking they should have been added to the quantum library.
 
 Our reason for not doing so was to make clear which code which was
 created within the week of the library generation and that code that
 needed to be added afterwards.
 
 */


import Foundation
import Quantum
public func spinWignerParity(theta: Double, phi: Double, spinSpace: StateSpace) -> MatrixOperator {
    let I = Complex(real: 0.0, imag: 1.0)
    assert (spinSpace.dimension == 2)
    
    let exp_i_simga_z_phi = cos(0.5 * phi) * spinSpace.identityOperator
    + (I * sin(0.5 * phi) ) * spinSpace.sigmaZ
    let exp_i_simga_y_theta = cos(0.5 * theta) * spinSpace.identityOperator
    + (I * sin(0.5 * theta) ) * spinSpace.sigmaY
    let spinParity = 0.5 * (spinSpace.identityOperator + sqrt(3.0) * spinSpace.sigmaZ )
    
    let displacementOperator = exp_i_simga_z_phi * exp_i_simga_y_theta
    let Pi = displacementOperator * spinParity * displacementOperator.hermitianAdjoint()
    return Pi
}

public func fieldWignerParity(q: Double, p: Double, fieldSpace: StateSpace) -> MatrixOperator  {
    let twoalpha = 2.0 * Complex(real: q, imag: p)
    let pre_factor = exp( twoalpha.modulus * twoalpha.modulus * -0.5 )
    // https://journals.aps.org/pra/abstract/10.1103/PhysRevA.50.4488 D(2a)Pi = D(a) Pi D^dag(a)
    let displacementOperator = Complex(real: pre_factor) *
    fieldSpace.exponentialOfScaledCreationOperator(scaleFactor: twoalpha) *
    fieldSpace.exponentialOfScaledAnnihilationOperator(scaleFactor: -(twoalpha.conjugate) )
    let Pi = displacementOperator * fieldSpace.parityOperator
    return Pi
}


// https://en.wikipedia.org/wiki/Box–Muller_transform
// Polar form
struct RandomGaussianZeroMeanUnitVariance {
    static var have_v = false
    static var radiusSquared = 0.0
    static var u = 0.0
    static var v = 0.0

    public static func nextReal() -> Double {

        if have_v {
            have_v = false
            return v * sqrt(-2.0 * log(radiusSquared) / radiusSquared )
        }

        repeat {
            u = Double.random(in: -1 ..< 1)
            v = Double.random(in: -1 ..< 1)
            radiusSquared = u * u + v * v
        } while radiusSquared >= 1.0

        have_v = true
        return u * sqrt(-2.0 * log(radiusSquared) / radiusSquared )
    }
    
    public static func nextComplex() -> Complex<Double> {
        let amplitude = Self.nextReal()
        let phase = Double.random(in: 0 ..< Double.pi)
        return Complex(real: amplitude * cos(phase), imag: amplitude * sin(phase) )
    }
}



extension MatrixOperator {
    public func plus_identity_scaled(by x: ComplexReal) -> MatrixOperator {
        
        var output = self
        
        for i in 0 ..< space.dimension {
            output[i, i] = output[i, i] + x
        }
        
        return output
    }
    public func plus_identity_scaled(by x: Real) -> MatrixOperator {
        
        var output = self
        
        for i in 0 ..< space.dimension {
            output[i, i] = output[i, i] + x
        }
        
        return output
    }
}

extension StateSpace {
    // https://en.wikipedia.org/wiki/Squeezed_coherent_state
    // |SMSV> in Single-mode squeezed states section
    public func makeSqueezedVacuumState(r: Double, phi: Double) -> StateVector {
        
        var output = Vector(in: self)
        
        typealias R = ScalarField.ScalarField
        
        let prefactor = R(1.0 / sqrt(cosh(r)))
        
        var factorialfactor = 1.0
        var squeezefactor = ScalarField(real: 1.0)

        output[0] = ScalarField(prefactor)
        
        var n = 1
        while 2 * n < self.dimension {
            let twon = 2.0 * Double(n)
            // avoid computing n! etc as this can get problematic
            factorialfactor = sqrt( twon * (twon - 1)  ) / twon
            squeezefactor   = ScalarField(modulus: -tanh(r), argument: phi)
            output[2 * n]   = output[2 * n - 2] * squeezefactor * ScalarField(real: factorialfactor)
            n = n + 1
        }
        
        return output
    }
}

//
//  Created by M J Everitt on 28/08/2022.
//
