document.querySelectorAll('a[href^="https://"]').forEach(el=>{
    el.setAttribute("title", "This link will exit this website.");
    el.setAttribute("target", "_blank");
    el.setAttribute("rel", "noopener noreferrer");
});
