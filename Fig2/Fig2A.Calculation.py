input_info = open("NewtonParameters_DropoutDisturbing.txt","r")
u_info = open("Parameters_u_DropoutDisturbing.txt","w")
d_info = open("Parameters_d_DropoutDisturbing.txt","w")
p_info = open("Parameters_p_DropoutDisturbing.txt","w")

count = 0

for each in input_info:
	each = each.strip().split("\t")
	if count==0:
		tmp_u = [each[1]]
		tmp_d = [each[2]]
		tmp_p = [each[3]]
	elif count%9==0:
		u_info.write("{}\t\n".format("\t".join(tmp_u)))
		d_info.write("{}\t\n".format("\t".join(tmp_d)))
		p_info.write("{}\t\n".format("\t".join(tmp_p)))
		tmp_u = [each[1]]
		tmp_d = [each[2]]
		tmp_p = [each[3]]
	else:
		tmp_u.append(each[1])
		tmp_d.append(each[2])
		tmp_p.append(each[3])

	count += 1

u_info.write("{}\t\n".format("\t".join(tmp_u)))
d_info.write("{}\t\n".format("\t".join(tmp_d)))
p_info.write("{}\t\n".format("\t".join(tmp_p)))

input_info.close()
u_info.close()
d_info.close()
p_info.close()