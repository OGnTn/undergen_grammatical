@tool
extends Node

class_name PoissonDiskSampler

var valid_points = 0
var invalid_points = 0

# Helper function to check if a point is valid
func is_valid_point(sampled_points, new_point: Vector3, radius: float) -> bool:
	for existing_point: Vector3 in sampled_points:
		var distance = existing_point.distance_to(new_point)
		if distance < radius:
			#print("Point is too close to existing point")
			invalid_points += 1
			return false
	valid_points += 1
	return true

func poisson_disk_sampling(points: Array, radius: float) -> Array:
	# Convert dictionary keys to an Array of Vector3
	var positions = points
	var sampled_points: Array = []
	
	# Iterate over points and add valid points to the sampled list
	for point in positions:
		if is_valid_point(sampled_points, point, radius):
			sampled_points.append(point)
	
	print("Valid points: ", valid_points)
	print("Invalid points: ", invalid_points)
	print("Sampled points: ", sampled_points.size())
	print("Total points: ", positions.size())
	return sampled_points
