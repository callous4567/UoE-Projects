using MPI
MPI.Init()

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
print("Rank is $rank")
print("\n")
print(MPI.universe_size())
print("\n")
N = MPI.Comm_size(comm) 
print(N)
print("\n")

recv_buf = Array{Float64}(undef, 2)
recv_req = MPI.Irecv!(recv_buf, mod(rank-1, N), 0,
                    comm)

send_buf = Float64[rank, rank]
send_req = MPI.Isend(send_buf, mod(rank+1, N), 0,
                    comm)


MPI.Waitall!([recv_req, send_req])
print("$rank: Received $recv_buf\n")