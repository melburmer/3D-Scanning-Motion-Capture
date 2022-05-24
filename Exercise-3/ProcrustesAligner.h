#pragma once
#include "SimpleMesh.h"

class ProcrustesAligner {
public:
	Matrix4f estimatePose(const std::vector<Vector3f>& sourcePoints, const std::vector<Vector3f>& targetPoints) {
		ASSERT(sourcePoints.size() == targetPoints.size() && "The number of source and target points should be the same, since every source point is matched with corresponding target point.");

		// We estimate the pose between source and target points using Procrustes algorithm.
		// Our shapes have the same scale, therefore we don't estimate scale. We estimated rotation and translation
		// from source points to target points.

		auto sourceMean = computeMean(sourcePoints);
		auto targetMean = computeMean(targetPoints);
		
		Matrix3f rotation = estimateRotation(sourcePoints, sourceMean, targetPoints, targetMean);
		Vector3f translation = computeTranslation(sourceMean, targetMean, rotation);

		// TODO: Compute the transformation matrix by using the computed rotation and translation.
		// You can access parts of the matrix with .block(start_row, start_col, num_rows, num_cols) = elements
		
		Matrix4f estimatedPose = Matrix4f::Identity();
		estimatedPose.block<3, 3>(0, 0) = rotation;
		estimatedPose.block<3, 1>(0, 3) = translation;
		return estimatedPose;
	}

private:
	Vector3f computeMean(const std::vector<Vector3f>& points) {
		// TODO: Compute the mean of input points.
		// Hint: You can use the .size() method to get the length of a vector.

		Vector3f mean = Vector3f::Zero();
		Vector3f total = Vector3f::Zero();

		for (unsigned i = 0; i < points.size(); ++i)
			total += points[i];

		mean = total / points.size();

		return mean;
	}

	Matrix3f estimateRotation(const std::vector<Vector3f>& sourcePoints, const Vector3f& sourceMean, const std::vector<Vector3f>& targetPoints, const Vector3f& targetMean) {
		// TODO: Estimate the rotation from source to target points, following the Procrustes algorithm.
		// To compute the singular value decomposition you can use JacobiSVD() from Eigen.
		// Hint: You can initialize an Eigen matrix with "MatrixXf m(num_rows,num_cols);" and access/modify parts of it using the .block() method (see above).

		Matrix3f rotation = Matrix3f::Identity();

		// compute zero centered data
		std::vector<Vector3f> sourcePoints_centered;
		for (unsigned i = 0; i < sourcePoints.size(); ++i)
			sourcePoints_centered.push_back(sourcePoints[i] - sourceMean);

		std::vector<Vector3f> targetPoints_centered;
		for (unsigned i = 0; i < targetPoints.size(); ++i)
			targetPoints_centered.push_back(targetPoints[i] - targetMean);

		// Calculate X.T*X.hat
		Matrix3f Matrix_Mul = Matrix3f::Zero();
		for (unsigned i = 0; i < sourcePoints.size(); ++i)
			Matrix_Mul += targetPoints[i] * sourcePoints[i].transpose();

		// SVD Decomposition
		Eigen::JacobiSVD<Matrix3f> svd(Matrix_Mul, Eigen::ComputeFullU | Eigen::ComputeFullV);
		// Calculate determinant of U*V.T
		double det{ (svd.matrixU() * svd.matrixV().transpose()).determinant() };
		if (det == -1)
		{
			// fix reflection
			Matrix3f a{
				{1,0,0},
				{0,1,0},
				{0,0,-1},
			};
			rotation = svd.matrixU() * a * svd.matrixV().transpose();
		}
		else
		{
			rotation = svd.matrixU() * svd.matrixV().transpose();
		}


        return rotation;
	}

	Vector3f computeTranslation(const Vector3f& sourceMean, const Vector3f& targetMean, const Matrix3f& rotation) {
		// TODO: Compute the translation vector from source to target points.

		Vector3f translation = Vector3f::Zero();
		translation = -1 * rotation * sourceMean + targetMean;
        return translation;
	}
};
