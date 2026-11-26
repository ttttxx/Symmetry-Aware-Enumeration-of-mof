import sys
def solve():
    data = sys.stdin.read().split()
    t = int(data[0])
    idx = 1
    results = []

    for _ in range(t):
        n = int(data[idx])
        m = int(data[idx+1])
        idx += 2

        prices = [int(data[idx+i]) for i in range(n)]
        idx += n

        no_discount = {0}
        use_discount = set()

        for a in prices:
            half = a//2
            new_use = set(use_discount)
            for s in no_discount :
                if s + half <= m:
                    new_use.add(s+half)
            use_discount =new_use

            new_no_discount = set(no_discount)
            for s in no_discount:
                if s + a <= m:
                    new_no_discount.add(s+a)
            no_discount = new_no_discount
        if m in no_discount or m in use_discount:
            results.append("YES")
        else:
            results.append("NO")
    print('\n'.join(results))

if __name__=="__main__":
    solve()