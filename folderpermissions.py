BASE_ROUTE = '/PHYLOISLAND'

def local(route: str) -> str:
    if BASE_ROUTE == '/':
        return route
    else:
        return join(BASE_ROUTE, route[1:])

def local_url_for(*args, **kwargs) -> str:
    new_url = local(url_for(*args, **kwargs))
    if new_url.count(BASE_ROUTE[1:]) == 1:
        fixed_url = new_url
        return new_url
    else:
        fixed_url = '/'.join(new_url.split('/')[2:])
    assert fixed_url.count(BASE_ROUTE[1:]) == 1, fixed_url
    return fixed_url

def local_redirect(*args, **kwargs) -> Any:
    return redirect(local_url_for(*args, **kwargs))