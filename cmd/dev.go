package main

import (
	"fmt"
	"github.com/zedseven/bch"
	"math/rand"
)

func corruptData(data *[]int, errors int) {
	f := randIndex(len(*data))
	for i := 0; i < errors; i++ {
		l, _ := f()
		(*data)[l] ^= 1
	}
}

func randIndex(count int) func() (int, error) {
	pool := makeRange(int64(count))
	return func() (int, error) {
		if count <= 0 {
			return -1, nil
		}

		j := rand.Int63n(int64(count)) //I'm aware this isn't crypto/rand, but I needed to be able to seed it

		count--

		p := int(pool[j])

		pool[j] = pool[count]
		pool = pool[:count]

		return p, nil
	}
}

func makeRange(max int64) []int64 {
	r := make([]int64, max)
	for i := range r {
		r[i] = int64(i)
	}
	return r
}

func main() {
	/*fmt.Println("Enter the code length: ")
	length := 0
	fmt.Scanf("%d", &length)
	m := uint8(math.Log2(float64(length))) + 1
	bch.readP(m)*/

	lol := []int{1, 1, 1, 0, 0, 1, 0, 0}
	fmt.Printf("Original data: %v\n", lol)
	data, config := bch.Encode(32, 4, &lol)
	fmt.Printf("Encoded data:   %v\n", data)
	corruptData(&data, 4)
	fmt.Printf("Corrupted data: %v\n", data)
	recv := bch.Decode(&data, config)
	fmt.Printf("Decoded data:   %v\n", recv)
}