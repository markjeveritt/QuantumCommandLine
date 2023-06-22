//
//  SparseSE.swift
//  QuantumCommandLine
//
//  Created by DOCTOR M J EVERITT on 25/08/2021.
//
//
//import Foundation
//import Quantum
//
//public class SparseSchrodingerEquationConstantHamiltoian: SimpleIntegrableSystem {
//
//
//    public typealias ScalarType = myChoiceOfNumberField
//    public typealias ComplexType = Complex<ScalarType>
//    public typealias IntegrableType = ColumnVector<ComplexType>
//
//    public var currentState: IntegrableType
//    public var currentIndependentVariable: ScalarType
//    public var time: ScalarType { return currentIndependentVariable }
//
//    let minusI = -ComplexType.I
//
//    public func computeDerivatives(x: ScalarType, y: IntegrableType) -> IntegrableType {
//        return ( minusI_times_HamiltomianOverHbar * y )
//    }
//
//    let minusI_times_HamiltomianOverHbar: SparseMatrixOperator<ComplexType>
//
//    public init(startTime: Double, initialstate: IntegrableType, hamiltonian: MatrixOperator<ComplexType>){
//        currentIndependentVariable = startTime
//        currentState = initialstate
//        minusI_times_HamiltomianOverHbar = SparseMatrixOperator(from: hamiltonian) * minusI
//    }
//    public static func times(vector: ColumnVector<ComplexType>, by: Double) -> ColumnVector<ComplexType> {
//        return vector * by
//    }
//
//}
